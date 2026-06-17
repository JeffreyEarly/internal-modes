classdef IMBasisSetExponentialStratification < IMInternalModesBasis
    % Evaluate exact basis sets for exponential stratification.
    %
    % `IMBasisSetExponentialStratification` stores exact rigid-bottom
    % `G`-formulation basis sets for rigid or free surfaces with
    % $$N^2(z)=N_0^2 e^{2z/b}$$ on domains with surface at $$z=0$$. Mode
    % roots are found with a local scanner and `fzero`, without requiring
    % an external spectral-function root finder at runtime.
    %
    % ```matlab
    % zDomain = [-5000 0];
    % N2 = @(z) (5.2e-3)^2*exp(2*z/1300);
    % evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
    % basisSet = IMBasisSetExponentialStratification(evp=evp, N0=5.2e-3, b=1300, zDomain=[-5000 0]);
    % G = basisSet.G(linspace(-5000,0,128).');
    % ```
    %
    % - Topic: Create basis sets
    % - Topic: Evaluate basis sets
    % - Topic: Inspect basis sets
    % - Topic: Developer topics
    % - Declaration: classdef IMBasisSetExponentialStratification < IMInternalModesBasis

    properties (SetAccess = private)
        % Surface buoyancy frequency $$N_0$$ in radians per second.
        %
        % - Topic: Inspect basis sets
        N0

        % Exponential e-folding depth $$b$$ in meters.
        %
        % - Topic: Inspect basis sets
        b

        % Dimensionless Bessel roots used to construct the retained modes.
        %
        % - Topic: Inspect basis sets
        roots

        % Phase speeds $$c_j=\sqrt{g h_j}$$ for each retained mode.
        %
        % - Topic: Inspect basis sets
        phaseSpeeds

        % Modal frequencies in radians per second.
        %
        % - Topic: Inspect basis sets
        frequencies

        % Sign applied to each raw mode so that $$F_j(0)>0$$.
        %
        % - Topic: Inspect basis sets
        signFactors

        % Internal analytical branch kind for each retained mode.
        %
        % Values are `"bessel"` for ordinary exponential modes and
        % `"null"` for the hydrostatic `F` null mode.
        %
        % - Topic: Inspect basis sets
        modeKinds
    end

    methods
        function self = IMBasisSetExponentialStratification(options)
            % Create an exact exponential-stratification basis set.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSetExponentialStratification(options)
            % - Parameter options.evp: supported rigid-endpoint `G` eigenvalue-problem descriptor
            % - Parameter options.N0: surface buoyancy frequency
            % - Parameter options.b: exponential e-folding depth
            % - Parameter options.zDomain: physical vertical domain with surface at zero
            % - Parameter options.nModes: number of retained modes
            % - Parameter options.normalization: active normalization rule; omitted uses the basis-set default
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: exact exponential-stratification basis set
            arguments
                options.evp = []
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.b (1,1) double {mustBePositive} = 1300
                options.zDomain (1,2) double = [-1 0]
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.normalization = []
                options.metadata struct = struct()
            end

            zDomain = sort(options.zDomain);
            IMBasisSetExponentialStratification.validateDomain(zDomain);
            N2 = @(z) options.N0*options.N0*exp(2*z/options.b);
            evp = IMBasisSetExponentialStratification.resolveEVP(options.evp, N2, zDomain);
            [f0, g] = IMBasisSetExponentialStratification.physicalConstants(evp);
            [h, roots, frequencies, modeNumber, modeKinds] = IMBasisSetExponentialStratification.solveSpectrum( ...
                evp, options.N0, options.b, zDomain, options.nModes, f0, g);
            phaseSpeeds = sqrt(g*h);
            signFactors = IMBasisSetExponentialStratification.surfaceSignFactors( ...
                options.N0, options.b, zDomain, frequencies, phaseSpeeds, g, modeKinds);
            eigenvalues = 1 ./ h;

            modeSelectionDiagnostics = struct();
            metadata = options.metadata;
            metadata.analyticalBasis = "exponentialStratification";

            self@IMInternalModesBasis(evp=evp, nativeModes=zeros(0,length(h)), ...
                eigenvalues=eigenvalues, h=h, modeNumber=modeNumber, modeSelectionDiagnostics=modeSelectionDiagnostics, ...
                normalization=options.normalization, metadata=metadata, zDomain=zDomain, ...
                N2=N2);
            self.N0 = options.N0;
            self.b = options.b;
            self.roots = roots;
            self.phaseSpeeds = phaseSpeeds;
            self.frequencies = frequencies;
            self.signFactors = signFactors;
            self.modeKinds = modeKinds;
        end
    end

    methods (Access = protected)
        function values = rawVariable(self, variable, z)
            variable = string(variable);
            if variable ~= "G" && variable ~= "F"
                self.unsupported("evaluate " + variable);
            end

            z = z(:);
            values = zeros(length(z), length(self.h));
            for iMode = 1:length(self.h)
                if self.modeKinds(iMode) == "null"
                    G = zeros(size(z));
                    F = ones(size(z));
                else
                    [G, F] = IMBasisSetExponentialStratification.rawModeValues( ...
                        z, self.N0, self.b, self.zDomain, self.frequencies(iMode), ...
                        self.phaseSpeeds(iMode), self.evp.g);
                end
                G = self.signFactors(iMode)*G;
                F = self.signFactors(iMode)*F;
                if variable == "G"
                    values(:,iMode) = G;
                else
                    values(:,iMode) = F;
                end
            end
        end

        function values = rawUz(self, z)
            z = z(:);
            values = zeros(length(z), length(self.h));
            for iMode = 1:length(self.h)
                if self.modeKinds(iMode) == "null"
                    values(:,iMode) = zeros(size(z));
                    continue;
                end
                [G, F] = IMBasisSetExponentialStratification.rawModeValues( ...
                    z, self.N0, self.b, self.zDomain, self.frequencies(iMode), ...
                    self.phaseSpeeds(iMode), self.evp.g);
                G = self.signFactors(iMode)*G;
                F = self.signFactors(iMode)*F;

                if self.evp.formulation == "G"
                    values(:,iMode) = F ./ self.h(iMode);
                else
                    values(:,iMode) = -(self.N2(z)/self.evp.g).*G;
                end
            end
        end
    end

    methods (Static, Access = private)
        function evp = resolveEVP(evp, N2, zDomain)
            if isempty(evp)
                evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
                return;
            end
            if ~isa(evp, "IMInternalModes")
                error("IMBasisSetExponentialStratification:InvalidEVP", ...
                    "The EVP must be an IMInternalModes instance.");
            end
            IMBasisSetExponentialStratification.validateEVPDomain(evp, zDomain);
        end

        function validateEVPDomain(evp, zDomain)
            tolerance = 100*eps(max([1 abs(evp.zDomain) abs(zDomain)]));
            if max(abs(evp.zDomain - zDomain)) > tolerance
                error("IMBasisSetExponentialStratification:DomainMismatch", ...
                    "The analytical basis zDomain must match evp.zDomain.");
            end
        end

        function validateDomain(zDomain)
            tolerance = 100*eps(max(1,max(abs(zDomain))));
            if abs(zDomain(2)) > tolerance
                error("IMBasisSetExponentialStratification:UnsupportedDomain", ...
                    "Exponential-stratification analytical modes currently require the surface at z=0.");
            end
        end

        function [f0, g] = physicalConstants(evp)
            f0 = evp.f0;
            g = evp.g;
            if ~(isscalar(f0) && isfinite(f0))
                error("IMBasisSetExponentialStratification:InvalidCoriolis", ...
                    "The Coriolis parameter must be finite.");
            end
            if ~(isscalar(g) && isfinite(g) && g > 0)
                error("IMBasisSetExponentialStratification:InvalidGravity", ...
                    "The gravitational acceleration must be positive.");
            end
        end

        function [h, roots, frequencies, modeNumber, modeKinds] = solveSpectrum(evp, N0, b, zDomain, nModes, f0, g)
            surfaceBoundary = IMBasisSetExponentialStratification.validateEVP(evp);
            evpName = string(evp.name);
            switch evpName
                case "waveModesAtWavenumber"
                    if ~isfield(evp.parameters, "k")
                        error("IMBasisSetExponentialStratification:UnsupportedEVP", ...
                            "A fixed-wavenumber EVP must include parameters.k.");
                    end
                    k = evp.parameters.k;
                    nInteriorModes = IMBasisSetExponentialStratification.nInteriorModes(nModes, surfaceBoundary);
                    roots = IMBasisSetExponentialStratification.rootsAtWavenumber( ...
                        k, N0, b, zDomain, nInteriorModes, f0, g, surfaceBoundary);
                    h = (b*N0./roots).^2/g;
                    if surfaceBoundary == "free"
                        root0 = IMBasisSetExponentialStratification.surfaceRootAtWavenumber(k, N0, b, zDomain, f0, g);
                        roots = [root0 roots];
                        h = [(b*N0/root0)^2/g h];
                    end
                    frequencies = sqrt(g*h*k*k + f0*f0);
                case "waveModesAtFrequency"
                    if ~isfield(evp.parameters, "omega")
                        error("IMBasisSetExponentialStratification:UnsupportedEVP", ...
                            "A fixed-frequency EVP must include parameters.omega.");
                    end
                    omega = evp.parameters.omega;
                    if omega >= N0 && surfaceBoundary == "free"
                        root0 = IMBasisSetExponentialStratification.surfaceRootAtFrequency(omega, N0, b, zDomain, g);
                        roots = root0;
                        h = (b*N0/root0)^2/g;
                    else
                        nInteriorModes = IMBasisSetExponentialStratification.nInteriorModes(nModes, surfaceBoundary);
                        [roots, eta] = IMBasisSetExponentialStratification.rootsAtFrequency( ...
                            omega, N0, b, zDomain, nInteriorModes, g, surfaceBoundary);
                        h = (b*eta./roots).^2/g;
                        if surfaceBoundary == "free"
                            root0 = IMBasisSetExponentialStratification.surfaceRootAtFrequency(omega, N0, b, zDomain, g);
                            roots = [root0 roots];
                            h = [(b*N0/root0)^2/g h];
                        end
                    end
                    frequencies = omega*ones(size(h));
                case "hydrostaticGModes"
                    omega = 0;
                    nInteriorModes = IMBasisSetExponentialStratification.nInteriorModes(nModes, surfaceBoundary);
                    [roots, eta] = IMBasisSetExponentialStratification.rootsAtFrequency( ...
                        omega, N0, b, zDomain, nInteriorModes, g, surfaceBoundary);
                    h = (b*eta./roots).^2/g;
                    if surfaceBoundary == "free"
                        root0 = IMBasisSetExponentialStratification.surfaceRootAtFrequency(omega, N0, b, zDomain, g);
                        roots = [root0 roots];
                        h = [(b*N0/root0)^2/g h];
                    end
                    frequencies = zeros(size(h));
                case "hydrostaticFModes"
                    omega = 0;
                    nInteriorModes = max(nModes - 1,0);
                    [rootsInterior, eta] = IMBasisSetExponentialStratification.rootsAtFrequency( ...
                        omega, N0, b, zDomain, nInteriorModes, g, "rigid");
                    hInterior = (b*eta./rootsInterior).^2/g;
                    roots = [NaN rootsInterior];
                    h = [Inf hInterior];
                    frequencies = zeros(size(h));
                otherwise
                    error("IMBasisSetExponentialStratification:UnsupportedEVP", ...
                        "Exponential stratification supports fixed-wavenumber, fixed-frequency, " + ...
                        "hydrostatic G, and hydrostatic F EVPs.");
            end
            if evpName == "hydrostaticFModes"
                modeKinds = ["null" repmat("bessel",1,length(h)-1)];
            else
                modeKinds = repmat("bessel",1,length(h));
            end
            [roots, h, frequencies, modeKinds] = IMBasisSetExponentialStratification.truncateModes( ...
                roots, h, frequencies, modeKinds, nModes);
            modeNumber = IMBasisSetExponentialStratification.modeNumberFor(surfaceBoundary, length(h));
            if evpName == "hydrostaticFModes"
                modeNumber = 0:(length(h)-1);
            end
        end

        function surfaceBoundary = validateEVP(evp)
            surfaceBoundary = IMBasisSetExponentialStratification.boundaryKind(evp.surfaceBoundary, evp.formulation);
            bottomBoundary = IMBasisSetExponentialStratification.boundaryKind(evp.bottomBoundary, evp.formulation);
            if evp.formulation == "F"
                if evp.name ~= "hydrostaticFModes"
                    error("IMBasisSetExponentialStratification:UnsupportedEVP", ...
                        "Exponential stratification supports F-formulation EVPs only for hydrostatic F modes.");
                end
                if surfaceBoundary ~= "rigid" || bottomBoundary ~= "rigid"
                    error("IMBasisSetExponentialStratification:UnsupportedBoundary", ...
                        "Hydrostatic F exponential-stratification modes currently require rigid surface and rigid bottom boundaries.");
                end
                return;
            end
            if evp.formulation ~= "G"
                error("IMBasisSetExponentialStratification:UnsupportedEVP", ...
                    "Exponential stratification supports only F or G formulations.");
            end
            if bottomBoundary ~= "rigid"
                error("IMBasisSetExponentialStratification:UnsupportedBoundary", ...
                    "Exponential stratification currently supports only rigid bottom boundaries.");
            end
            if ~ismember(surfaceBoundary, ["rigid", "free"])
                error("IMBasisSetExponentialStratification:UnsupportedBoundary", ...
                    "Exponential stratification currently supports rigid or free surface boundaries.");
            end
        end

        function kind = boundaryKind(boundary, formulation)
            formulation = string(formulation);
            if formulation == "G"
                if IMBasisSetExponentialStratification.isCondition(boundary, [1 0 0 0])
                    kind = "rigid";
                elseif IMBasisSetExponentialStratification.isCondition(boundary, [0 1 1 0])
                    kind = "free";
                elseif IMBasisSetExponentialStratification.isCondition(boundary, [0 1 0 0])
                    kind = "noSlip";
                else
                    kind = "custom";
                end
            else
                if IMBasisSetExponentialStratification.isCondition(boundary, [0 1 0 0])
                    kind = "rigid";
                elseif IMBasisSetExponentialStratification.isCondition(boundary, [1 0 0 0])
                    kind = "noSlip";
                else
                    kind = "custom";
                end
            end
        end

        function tf = isCondition(boundary, coefficients)
            actual = [boundary.a boundary.b boundary.c boundary.d];
            tolerance = 100*eps(max(1,max(abs([actual coefficients]))));
            tf = max(abs(actual - coefficients)) <= tolerance;
        end

        function roots = rootsAtWavenumber(k, N0, b, zDomain, nModes, f0, g, surfaceBoundary)
            if nModes == 0
                roots = zeros(1,0);
                return;
            end
            D = diff(zDomain);
            epsilon = f0/N0;
            lambda = k*b;
            x_lf = @(j,lambdaValue) (j - 1/4)*pi + lambdaValue*pi/2;
            x_hf = @(j,lambdaValue) lambdaValue.*(1 + 0.5*(3*pi*(4*j - 1)./(lambdaValue*8*sqrt(2))).^(2/3));
            if lambda < 2*(1 - 1/4)*1e-1
                lowerBound = x_lf(1, lambda);
            else
                lowerBound = x_hf(1, lambda);
            end
            if lambda < (IMBasisSetExponentialStratification.nInitialSearchModes() - 1/4)
                upperBound = x_lf(1.1*IMBasisSetExponentialStratification.nInitialSearchModes(), lambda);
            else
                upperBound = x_hf(5*IMBasisSetExponentialStratification.nInitialSearchModes(), lambda);
            end

            nu = @(x) sqrt(epsilon*epsilon*x.*x + lambda*lambda);
            s = @(x) x;
            roots = IMBasisSetExponentialStratification.findEnoughRoots(nu, s, ...
                [lowerBound upperBound], exp(-D/b), nModes, surfaceBoundary, N0, b, g);
        end

        function [roots, eta] = rootsAtFrequency(omega, N0, b, zDomain, nModes, g, surfaceBoundary)
            D = diff(zDomain);
            expMinusDOverB = exp(-D/b);
            if omega >= N0
                error("IMBasisSetExponentialStratification:UnsupportedFrequency", ...
                    "Fixed-frequency exponential modes require omega < N0.");
            end
            if omega > N0*expMinusDOverB
                eta = (sqrt(N0*N0 - omega*omega) - omega*acos(omega/N0))/pi;
            else
                eta = (sqrt(N0*N0 - omega*omega) ...
                    - sqrt(N0*N0*exp(-2*D/b) - omega*omega) ...
                    - omega*acos(omega/N0) + omega*acos(omega*exp(D/b)/N0))/pi;
            end

            nu = @(x) omega*x/eta;
            s = @(x) N0*x/eta;
            if nModes == 0
                roots = zeros(1,0);
                return;
            end
            roots = IMBasisSetExponentialStratification.findEnoughRoots(nu, s, ...
                [0.5 nModes+1], expMinusDOverB, nModes, surfaceBoundary, N0, b, g);
        end

        function root = surfaceRootAtWavenumber(k, N0, b, zDomain, f0, g)
            D = diff(zDomain);
            kSafe = max(k,1e-15);
            rootEstimate = b*N0/sqrt(g*tanh(kSafe*D)/kSafe);
            epsilon = f0/N0;
            lambda = k*b;
            nu = @(x) sqrt(epsilon*epsilon*x.*x + lambda*lambda);
            s = @(x) x;
            root = IMBasisSetExponentialStratification.findFirstRoot( ...
                nu, s, [0.95 1.05]*rootEstimate, exp(-D/b), "free", N0, b, g);
        end

        function root = surfaceRootAtFrequency(omega, N0, b, zDomain, g)
            D = diff(zDomain);
            rootEstimate = max(b*N0*omega/g, b*N0/sqrt(g*D));
            nu = @(x) omega*x/N0;
            s = @(x) x;
            root = IMBasisSetExponentialStratification.findFirstRoot( ...
                nu, s, [0.95 1.5]*rootEstimate, exp(-D/b), "free", N0, b, g);
        end

        function root = findFirstRoot(nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, g)
            roots = IMBasisSetExponentialStratification.findRootsInRange( ...
                nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, g);
            if isempty(roots)
                error("IMBasisSetExponentialStratification:RootSearchFailed", ...
                    "Could not find the exponential free-surface root.");
            end
            root = roots(1);
        end

        function roots = findEnoughRoots(nu, s, bounds, expMinusDOverB, nModes, surfaceBoundary, N0, b, g)
            roots = IMBasisSetExponentialStratification.findRootsInRange( ...
                nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, g);
            iteration = 0;
            while length(roots) < nModes && iteration < 20
                iteration = iteration + 1;
                if length(roots) >= 2
                    dx = roots(end) - roots(end-1);
                else
                    dx = diff(bounds)/IMBasisSetExponentialStratification.nInitialSearchModes();
                end
                dx = max(dx, eps(max(1,bounds(2))));
                oldUpperBound = bounds(2);
                nSearchModes = min(IMBasisSetExponentialStratification.nInitialSearchModes(), ...
                    nModes - length(roots) + 1);
                bounds = [oldUpperBound oldUpperBound + dx*nSearchModes];
                roots = IMBasisSetExponentialStratification.deduplicateRoots( ...
                    [roots; IMBasisSetExponentialStratification.findRootsInRange( ...
                    nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, g)]);
            end
            if length(roots) < nModes
                error("IMBasisSetExponentialStratification:RootSearchFailed", ...
                    "Could not find %d exponential-stratification roots.", nModes);
            end
            roots = reshape(roots(1:nModes),1,[]);
        end

        function roots = findRootsInRange(nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, g)
            x = linspace(bounds(1), bounds(2), IMBasisSetExponentialStratification.nInitialSearchModes()).';
            intervals = zeros(0,2);

            bigNuIndices = find(nu(x) >= s(x)*expMinusDOverB);
            if ~isempty(bigNuIndices) && length(bigNuIndices) > 1
                intervals(end+1,:) = [x(min(bigNuIndices)) x(max(bigNuIndices))];
            end

            smallNuIndices = find(nu(x) < s(x)*expMinusDOverB);
            if ~isempty(smallNuIndices) && length(smallNuIndices) > 1
                if ~isempty(bigNuIndices)
                    intervals(end+1,:) = [x(max(bigNuIndices)) x(max(smallNuIndices))];
                else
                    intervals(end+1,:) = [x(min(smallNuIndices)) x(max(smallNuIndices))];
                end
            end

            roots = zeros(0,1);
            for iInterval = 1:size(intervals,1)
                interval = intervals(iInterval,:);
                if interval(2) <= interval(1)
                    continue;
                end
                if iInterval == 1 && ~isempty(bigNuIndices)
                    residual = @(x) IMBasisSetExponentialStratification.bigNuResidual( ...
                        nu, s, expMinusDOverB, x, surfaceBoundary, N0, b, g);
                else
                    residual = @(x) IMBasisSetExponentialStratification.smallNuResidual( ...
                        nu, s, expMinusDOverB, x, surfaceBoundary, N0, b, g);
                end
                roots = [roots; IMBasisSetExponentialStratification.scanRoots(residual, interval)];
            end
            roots = IMBasisSetExponentialStratification.deduplicateRoots(roots);
        end

        function values = smallNuResidual(nu, s, expMinusDOverB, x, surfaceBoundary, N0, b, g)
            [A, B] = IMBasisSetExponentialStratification.surfaceCoefficients(nu, s, x, surfaceBoundary, N0, b, g);
            values = A.*besselj(nu(x), expMinusDOverB*s(x)) ...
                + B.*bessely(nu(x), expMinusDOverB*s(x));
        end

        function values = bigNuResidual(nu, s, expMinusDOverB, x, surfaceBoundary, N0, b, g)
            [A, B] = IMBasisSetExponentialStratification.surfaceCoefficients(nu, s, x, surfaceBoundary, N0, b, g);
            denominator = bessely(nu(x), expMinusDOverB*s(x));
            values = (A./denominator).*besselj(nu(x), expMinusDOverB*s(x)) + B;
        end

        function [A, B] = surfaceCoefficients(nu, s, x, surfaceBoundary, N0, b, g)
            order = nu(x);
            argument = s(x);
            switch surfaceBoundary
                case "rigid"
                    A = bessely(order, argument);
                    B = -besselj(order, argument);
                case "free"
                    surfaceCoefficient = b*N0*N0/(2*g);
                    A = bessely(order, argument) - (surfaceCoefficient./argument) ...
                        .*(bessely(order-1, argument) - bessely(order+1, argument));
                    B = -besselj(order, argument) + (surfaceCoefficient./argument) ...
                        .*(besselj(order-1, argument) - besselj(order+1, argument));
                otherwise
                    error("IMBasisSetExponentialStratification:UnsupportedBoundary", ...
                        "Unsupported exponential surface boundary ""%s"".", surfaceBoundary);
            end
        end

        function roots = scanRoots(residual, interval)
            nSamples = max(2048, ceil(64*diff(interval)));
            x = linspace(interval(1), interval(2), nSamples).';
            y = IMBasisSetExponentialStratification.safeEvaluate(residual, x);
            roots = zeros(0,1);
            finite = isfinite(y);
            branchScale = median(abs(y(finite)));
            if isempty(branchScale) || ~isfinite(branchScale) || branchScale == 0
                branchScale = 1;
            end
            residualTolerance = 1e-6*max(1,branchScale);

            for iSample = 1:(length(x)-1)
                if ~finite(iSample) || ~finite(iSample+1)
                    continue;
                end
                if y(iSample) == 0
                    roots(end+1,1) = x(iSample);
                elseif sign(y(iSample)) == sign(y(iSample+1))
                    continue;
                else
                    try
                        root = fzero(residual, [x(iSample) x(iSample+1)]);
                        residualAtRoot = IMBasisSetExponentialStratification.safeEvaluate(residual, root);
                        if isfinite(root) && isfinite(residualAtRoot) && abs(residualAtRoot) <= residualTolerance
                            roots(end+1,1) = root;
                        end
                    catch
                    end
                end
            end
            roots = IMBasisSetExponentialStratification.deduplicateRoots(roots);
        end

        function values = safeEvaluate(residual, x)
            try
                rawValues = residual(x);
            catch
                rawValues = arrayfun(residual, x);
            end
            values = real(rawValues);
            values(abs(imag(rawValues)) > 1e-10*max(1,abs(values))) = NaN;
        end

        function roots = deduplicateRoots(roots)
            roots = sort(real(roots(:)));
            roots = roots(isfinite(roots));
            if isempty(roots)
                return;
            end
            tolerance = 1e-8*max(1,max(abs(roots)));
            keep = [true; abs(diff(roots)) > tolerance];
            roots = roots(keep);
        end

        function [G, F] = rawModeValues(z, N0, b, zDomain, omega, c, g)
            D = diff(zDomain);
            nu = b*omega/c;
            argument = (b*N0/c)*exp(z(:)/b);
            bottomArgument = (b*N0/c)*exp(-D/b);
            alpha = besselj(nu, bottomArgument)/bessely(nu, bottomArgument);
            if abs(alpha) < 1e-15
                G = besselj(nu, argument);
                F = (N0*exp(z(:)/b)*c/(2*g)).*(besselj(nu-1, argument) - besselj(nu+1, argument));
            else
                G = besselj(nu, argument) - alpha*bessely(nu, argument);
                F = (N0*exp(z(:)/b)*c/(2*g)) ...
                    .*((besselj(nu-1, argument) - besselj(nu+1, argument)) ...
                    - alpha*(bessely(nu-1, argument) - bessely(nu+1, argument)));
            end
        end

        function signFactors = surfaceSignFactors(N0, b, zDomain, frequencies, phaseSpeeds, g, modeKinds)
            signFactors = ones(1,length(frequencies));
            for iMode = 1:length(frequencies)
                if modeKinds(iMode) == "null"
                    continue;
                end
                [~, FSurface] = IMBasisSetExponentialStratification.rawModeValues( ...
                    0, N0, b, zDomain, frequencies(iMode), phaseSpeeds(iMode), g);
                if FSurface < 0
                    signFactors(iMode) = -1;
                end
            end
        end

        function value = nInitialSearchModes()
            value = 128;
        end

        function nModes = nInteriorModes(totalModes, surfaceBoundary)
            if surfaceBoundary == "free"
                nModes = max(0,totalModes - 1);
            else
                nModes = totalModes;
            end
        end

        function [roots, h, frequencies, modeKinds] = truncateModes(roots, h, frequencies, modeKinds, nModes)
            keep = 1:min(nModes,length(h));
            roots = roots(keep);
            h = h(keep);
            frequencies = frequencies(keep);
            modeKinds = modeKinds(keep);
        end

        function modeNumber = modeNumberFor(surfaceBoundary, nModes)
            if surfaceBoundary == "free"
                modeNumber = [-1 1:(nModes-1)];
            else
                modeNumber = 1:nModes;
            end
        end
    end
end
