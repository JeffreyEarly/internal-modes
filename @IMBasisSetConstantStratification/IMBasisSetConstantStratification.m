classdef IMBasisSetConstantStratification < IMBasisSet
    % Evaluate exact v2 basis sets for constant stratification.
    %
    % `IMBasisSetConstantStratification` stores exact
    % constant-stratification `G` basis sets for $$N^2(z)=N_0^2$$. The class
    % implements the same `evaluate` and normalization contract as numerical
    % v2 basis sets, without storing a solver reference.
    %
    % ```matlab
    % evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4);
    % basisSet = IMBasisSetConstantStratification(evp=evp, N0=5.2e-3, zDomain=[-5000 0]);
    % G = basisSet.evaluate("G", linspace(-5000,0,128).');
    % ```
    %
    % - Topic: Create basis sets
    % - Topic: Evaluate basis sets
    % - Topic: Inspect basis sets
    % - Topic: Developer topics
    % - Declaration: classdef IMBasisSetConstantStratification < IMBasisSet

    properties (SetAccess = private)
        % Constant buoyancy frequency $$N_0$$ in radians per second.
        %
        % - Topic: Inspect basis sets
        N0

        % Vertical wavenumbers for each retained mode.
        %
        % - Topic: Inspect basis sets
        verticalWavenumbers

        % Analytical branch type for each retained mode.
        %
        % Values are `"baroclinic"`, `"linear"`, `"hyperbolic"`, or
        % `"trig"`.
        %
        % - Topic: Inspect basis sets
        solutionTypes

        % True for the free-surface barotropic branch.
        %
        % - Topic: Inspect basis sets
        isBarotropic
    end

    properties (Access = private)
        % Baroclinic mode number associated with each retained column.
        baroclinicNumbers
    end

    methods
        function self = IMBasisSetConstantStratification(options)
            % Create an exact constant-stratification basis set.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSetConstantStratification(options)
            % - Parameter options.evp: wave-mode or hydrostatic `G` eigenvalue-problem descriptor
            % - Parameter options.N0: constant buoyancy frequency
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.nModes: number of retained modes
            % - Parameter options.normalization: active normalization; omitted uses the EVP default
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: exact constant-stratification basis set
            arguments
                options.evp IMEigenvalueProblem
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.zDomain (1,2) double = [-1 0]
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.normalization = []
                options.metadata struct = struct()
            end

            zDomain = sort(options.zDomain);
            [f0, g] = IMBasisSetConstantStratification.physicalConstants(options.evp);
            [h, verticalWavenumbers, solutionTypes, isBarotropic, baroclinicNumbers, modeIndex] = ...
                IMBasisSetConstantStratification.solveSpectrum(options.evp, options.N0, zDomain, options.nModes, f0, g);
            eigenvalues = 1 ./ h;
            context.N2 = @(z) options.N0*options.N0*ones(size(z));
            context.dzLogN2 = @(z) zeros(size(z));
            context.f0 = f0;
            context.g = g;
            context.zDomain = zDomain;
            context.coordinateKind = "constantStratification";
            index = options.evp.indexPolicy.classify(eigenvalues(:), context);
            metadata = options.metadata;
            metadata.analyticalBasis = "constantStratification";

            self@IMBasisSet(evp=options.evp, nativeModes=zeros(0,length(h)), ...
                eigenvalues=eigenvalues, h=h, modeIndex=modeIndex, index=index, normalization=options.normalization, ...
                metadata=metadata, zDomain=zDomain, N2Function=@(z) options.N0*options.N0*ones(size(z)));
            self.N0 = options.N0;
            self.verticalWavenumbers = verticalWavenumbers;
            self.solutionTypes = solutionTypes;
            self.isBarotropic = isBarotropic;
            self.baroclinicNumbers = baroclinicNumbers;
        end

        function factors = normalizationFactors(self, normalization)
            % Return EVP-defined normalization factors.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: factors = normalizationFactors(basisSet,normalization)
            % - Parameter normalization: normalization convention
            % - Returns factors: row vector of normalization factors
            arguments
                self IMBasisSetConstantStratification
                normalization = self.normalization
            end

            factors = normalizationFactors@IMBasisSet(self, normalization);
        end
    end

    methods (Hidden)
        function factor = weightedNormFactor(self, component, iMode, innerWeight, surfaceWeight, bottomWeight)
            % Return an exact constant-stratification norm factor.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = weightedNormFactor(basisSet,component,iMode,innerWeight,surfaceWeight,bottomWeight)
            % - Parameter component: component name
            % - Parameter iMode: mode index
            % - Parameter innerWeight: interior weight coefficient
            % - Parameter surfaceWeight: surface endpoint weight
            % - Parameter bottomWeight: bottom endpoint weight
            % - Returns factor: divisor from the requested quadratic form
            component = string(component);
            if component ~= "G" && component ~= "F"
                factor = weightedNormFactor@IMBasisSet(self, component, iMode, innerWeight, surfaceWeight, bottomWeight);
                return;
            end

            context = self.context();
            sampleZ = linspace(self.zDomain(1), self.zDomain(2), 3).';
            weightSamples = IMOperator.evaluateCoefficient(innerWeight, sampleZ, context);
            weightTolerance = 100*eps(max(1,max(abs(weightSamples))));
            if any(~isfinite(weightSamples)) || max(abs(weightSamples - weightSamples(1))) > weightTolerance
                factor = weightedNormFactor@IMBasisSet(self, component, iMode, innerWeight, surfaceWeight, bottomWeight);
                return;
            end
            weight = weightSamples(1);

            [interiorIntegral, surfaceValue, bottomValue] = self.componentQuadraticPieces(component, iMode);
            normValue = weight*interiorIntegral + surfaceWeight*surfaceValue*surfaceValue + bottomWeight*bottomValue*bottomValue;
            factor = sqrt(abs(normValue));
        end

        function factor = maxAbsFactor(self, component, iMode)
            % Return an exact constant-stratification maximum-amplitude factor.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = maxAbsFactor(basisSet,component,iMode)
            % - Parameter component: component name
            % - Parameter iMode: mode index
            % - Returns factor: maximum absolute component amplitude
            component = string(component);
            if component ~= "G" && component ~= "F"
                factor = maxAbsFactor@IMBasisSet(self, component, iMode);
                return;
            end

            Lz = diff(self.zDomain);
            k_z = self.verticalWavenumbers(iMode);
            hMode = self.h(iMode);
            if ~self.isBarotropic(iMode)
                if component == "G"
                    factor = 1;
                else
                    factor = abs(hMode*k_z);
                end
                return;
            end

            switch string(self.solutionTypes(iMode))
                case "linear"
                    factor = Lz;
                case "hyperbolic"
                    if component == "G"
                        factor = abs(sinh(k_z*Lz));
                    else
                        factor = abs(hMode*k_z*cosh(k_z*Lz));
                    end
                case "trig"
                    if component == "G"
                        factor = abs(sin(k_z*Lz));
                    else
                        factor = abs(hMode*k_z);
                    end
                otherwise
                    factor = maxAbsFactor@IMBasisSet(self, component, iMode);
            end
        end
    end

    methods (Access = protected)
        function values = rawComponent(self, component, z)
            component = string(component);
            if component ~= "G" && component ~= "F"
                self.unsupported("evaluate " + component);
            end

            z = z(:);
            s = z - self.zDomain(1);
            values = zeros(length(z), length(self.h));
            for iMode = 1:length(self.h)
                k_z = self.verticalWavenumbers(iMode);
                hMode = self.h(iMode);
                if self.isBarotropic(iMode)
                    [G, F] = self.rawBarotropicMode(self.solutionTypes(iMode), k_z, hMode, s);
                else
                    signValue = (-1)^self.baroclinicNumbers(iMode);
                    G = signValue*sin(k_z*s);
                    F = signValue*hMode*k_z*cos(k_z*s);
                end

                if component == "G"
                    values(:,iMode) = G;
                else
                    values(:,iMode) = F;
                end
            end
        end
    end

    methods (Access = private)
        function [interiorIntegral, surfaceValue, bottomValue] = componentQuadraticPieces(self, component, iMode)
            Lz = diff(self.zDomain);
            k_z = self.verticalWavenumbers(iMode);
            hMode = self.h(iMode);
            if ~self.isBarotropic(iMode)
                [interiorIntegral, surfaceValue, bottomValue] = ...
                    self.trigonometricQuadraticPieces(component, k_z, hMode, Lz);
                return;
            end

            switch string(self.solutionTypes(iMode))
                case "linear"
                    if component == "G"
                        interiorIntegral = Lz*Lz*Lz/3;
                        surfaceValue = Lz;
                        bottomValue = 0;
                    else
                        interiorIntegral = Lz*Lz*Lz;
                        surfaceValue = Lz;
                        bottomValue = Lz;
                    end
                case "hyperbolic"
                    if component == "G"
                        interiorIntegral = sinh(2*k_z*Lz)/(4*k_z) - Lz/2;
                        surfaceValue = sinh(k_z*Lz);
                        bottomValue = 0;
                    else
                        interiorIntegral = hMode*hMode*k_z*k_z*(Lz/2 + sinh(2*k_z*Lz)/(4*k_z));
                        surfaceValue = hMode*k_z*cosh(k_z*Lz);
                        bottomValue = hMode*k_z;
                    end
                case "trig"
                    [interiorIntegral, surfaceValue, bottomValue] = ...
                        self.trigonometricQuadraticPieces(component, k_z, hMode, Lz);
                otherwise
                    error("IMBasisSetConstantStratification:InvalidSolutionType", ...
                        "Unknown solution type ""%s"".", string(self.solutionTypes(iMode)));
            end
        end

        function [interiorIntegral, surfaceValue, bottomValue] = trigonometricQuadraticPieces(~, component, k_z, hMode, Lz)
            if abs(k_z) <= eps
                interiorIntegral = 0;
                surfaceValue = 0;
                bottomValue = 0;
                return;
            end

            if component == "G"
                interiorIntegral = Lz/2 - sin(2*k_z*Lz)/(4*k_z);
                surfaceValue = sin(k_z*Lz);
                bottomValue = 0;
            else
                interiorIntegral = hMode*hMode*k_z*k_z*(Lz/2 + sin(2*k_z*Lz)/(4*k_z));
                surfaceValue = hMode*k_z*cos(k_z*Lz);
                bottomValue = hMode*k_z;
            end
        end

        function [G, F] = rawBarotropicMode(self, solutionType, k_z, hMode, s)
            Lz = diff(self.zDomain);
            switch string(solutionType)
                case "linear"
                    G = s;
                    F = Lz*ones(size(s));
                case "hyperbolic"
                    G = sinh(k_z*s);
                    F = hMode*k_z*cosh(k_z*s);
                case "trig"
                    G = sin(k_z*s);
                    F = hMode*k_z*cos(k_z*s);
                otherwise
                    error("IMBasisSetConstantStratification:InvalidSolutionType", ...
                        "Unknown solution type ""%s"".", string(solutionType));
            end
        end
    end

    methods (Static, Access = private)
        function [f0, g] = physicalConstants(evp)
            f0 = evp.f0;
            g = evp.g;
            if ~(isscalar(f0) && isfinite(f0))
                error("IMBasisSetConstantStratification:InvalidCoriolis", ...
                    "The Coriolis parameter must be finite.");
            end
            if ~(isscalar(g) && isfinite(g) && g > 0)
                error("IMBasisSetConstantStratification:InvalidGravity", ...
                    "The gravitational acceleration must be positive.");
            end
        end

        function [h, verticalWavenumbers, solutionTypes, isBarotropic, baroclinicNumbers, modeIndex] = solveSpectrum(evp, N0, zDomain, nModes, f0, g)
            IMBasisSetConstantStratification.validateEVP(evp);
            D = diff(zDomain);
            problemType = string(evp.parameters.problemType);
            upperBoundary = string(evp.parameters.upperBoundary);

            switch problemType
                case "waveModesAtWavenumber"
                    [hBaroclinic, k_zBaroclinic] = IMBasisSetConstantStratification.baroclinicAtWavenumber(evp.parameters.k, N0, D, nModes, f0, g, upperBoundary);
                    if upperBoundary == "free"
                        [h0, k_z0, solutionType0] = IMBasisSetConstantStratification.barotropicAtWavenumber(evp.parameters.k, N0, D, f0, g);
                        h = [h0 hBaroclinic(1:end-1)];
                        verticalWavenumbers = [k_z0 k_zBaroclinic(1:end-1)];
                        solutionTypes = IMBasisSetConstantStratification.freeSolutionTypes(solutionType0, nModes);
                        isBarotropic = [true false(1,nModes-1)];
                        baroclinicNumbers = [0 1:(nModes-1)];
                    else
                        h = hBaroclinic;
                        verticalWavenumbers = k_zBaroclinic;
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBarotropic = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "waveModesAtFrequency"
                    [hBaroclinic, k_zBaroclinic] = IMBasisSetConstantStratification.baroclinicAtFrequency(evp.parameters.omega, N0, D, nModes, upperBoundary, g);
                    if upperBoundary == "free"
                        [h0, k_z0, solutionType0] = IMBasisSetConstantStratification.barotropicAtFrequency(evp.parameters.omega, N0, D, g);
                        h = [h0 hBaroclinic(1:end-1)];
                        verticalWavenumbers = [k_z0 k_zBaroclinic(1:end-1)];
                        solutionTypes = IMBasisSetConstantStratification.freeSolutionTypes(solutionType0, nModes);
                        isBarotropic = [true false(1,nModes-1)];
                        baroclinicNumbers = [0 1:(nModes-1)];
                    else
                        h = hBaroclinic;
                        verticalWavenumbers = k_zBaroclinic;
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBarotropic = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "hydrostaticGModes"
                    [h, verticalWavenumbers] = IMBasisSetConstantStratification.baroclinicAtFrequency(0, N0, D, nModes, upperBoundary, g);
                    solutionTypes = repmat("baroclinic",1,nModes);
                    isBarotropic = false(1,nModes);
                    baroclinicNumbers = 1:nModes;
                otherwise
                    error("IMBasisSetConstantStratification:UnsupportedEVP", ...
                        "Constant stratification supports fixed-wavenumber, fixed-frequency, and hydrostatic G EVPs.");
            end
            modeIndex = baroclinicNumbers;
        end

        function [h, k_z] = baroclinicAtWavenumber(k, N0, D, nModes, f0, g, upperBoundary)
            k_z = (1:nModes)*pi/D;
            if upperBoundary == "free"
                for iMode = 1:nModes
                    f = @(xi) (xi + iMode*pi)*(N0*N0 - f0*f0)*D - g*(k*k*D*D + (xi + iMode*pi)*(xi + iMode*pi))*tan(xi);
                    k_z(iMode) = k_z(iMode) + fzero(f, 0)/D;
                end
            end
            h = (N0*N0 - f0*f0)./(g*(k*k + k_z.*k_z));
        end

        function solutionTypes = freeSolutionTypes(barotropicType, nModes)
            solutionTypes = repmat("baroclinic",1,nModes);
            solutionTypes(1) = barotropicType;
        end

        function [h, k_z] = baroclinicAtFrequency(omega, N0, D, nModes, upperBoundary, g)
            k_z = (1:nModes)*pi/D;
            if upperBoundary == "free"
                for iMode = 1:nModes
                    f = @(xi) g*tan(xi)/(xi + iMode*pi) - (N0*N0 - omega*omega)*D/((xi + iMode*pi)^2);
                    k_z(iMode) = k_z(iMode) + fzero(f, 0)/D;
                end
            end
            h = (N0*N0 - omega*omega)./(g*k_z.*k_z);
            if omega >= N0
                k_z = zeros(size(k_z));
                h = zeros(size(h));
            end
        end

        function [h0, k_z, solutionType] = barotropicAtWavenumber(k, N0, D, f0, g)
            kStar = sqrt((N0*N0 - f0*f0)/(g*D));
            if abs(k - kStar)/kStar < 1e-6
                solutionType = "linear";
                h0 = D;
                k_z = 0;
            elseif k > kStar
                solutionType = "hyperbolic";
                f = @(q) D*(N0*N0 - f0*f0) - (1./q).*(g*(k*k*D*D - q.*q)).*tanh(q);
                kInitial = sqrt(k*k*D*D - D*(N0*N0 - f0*f0)/g);
                k_z = fzero(f, kInitial)/D;
                h0 = (N0*N0 - f0*f0)/(g*(k*k - k_z*k_z));
            else
                solutionType = "trig";
                f = @(q) D*(N0*N0 - f0*f0) - (1./q).*(g*(k*k*D*D + q.*q)).*tan(q);
                kInitial = sqrt(-k*k*D*D + D*(N0*N0 - f0*f0)/g);
                k_z = fzero(f, kInitial)/D;
                h0 = (N0*N0 - f0*f0)/(g*(k*k + k_z*k_z));
            end
        end

        function [h0, k_z, solutionType] = barotropicAtFrequency(omega, N0, D, g)
            if abs(omega - N0)/N0 < 1e-6
                solutionType = "linear";
                h0 = D;
                k_z = 0;
            elseif omega > N0
                solutionType = "hyperbolic";
                f = @(q) D*(omega*omega - N0*N0) - g*q.*tanh(q);
                k_z = fzero(f, sqrt(D*(omega*omega - N0*N0)/g))/D;
                h0 = (omega*omega - N0*N0)/(g*k_z*k_z);
            else
                solutionType = "trig";
                f = @(q) D*(N0*N0 - omega*omega) - g*q.*tan(q);
                k_z = fzero(f, sqrt(D*(N0*N0 - omega*omega)/g))/D;
                h0 = (N0*N0 - omega*omega)/(g*k_z*k_z);
            end
        end

        function validateEVP(evp)
            if ~isfield(evp.parameters, "problemType")
                error("IMBasisSetConstantStratification:UnsupportedEVP", ...
                    "The EVP must include v2 problem parameters.");
            end
            if ~isfield(evp.parameters, "upperBoundary") || ~isfield(evp.parameters, "lowerBoundary")
                error("IMBasisSetConstantStratification:UnsupportedEVP", ...
                    "The EVP must include upper and lower boundary metadata.");
            end
            lowerBoundary = string(evp.parameters.lowerBoundary);
            upperBoundary = string(evp.parameters.upperBoundary);
            if ~ismember(lowerBoundary, ["rigid", "rigidLid", "freeSlip"])
                error("IMBasisSetConstantStratification:UnsupportedBoundary", ...
                    "Unsupported constant-stratification lower boundary ""%s"".", lowerBoundary);
            end
            if ~ismember(upperBoundary, ["rigid", "rigidLid", "free"])
                error("IMBasisSetConstantStratification:UnsupportedBoundary", ...
                    "Unsupported constant-stratification upper boundary ""%s"".", upperBoundary);
            end
        end
    end
end
