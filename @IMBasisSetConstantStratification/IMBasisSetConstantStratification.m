classdef IMBasisSetConstantStratification < IMInternalModesBasis
    % Evaluate exact basis sets for constant stratification.
    %
    % `IMBasisSetConstantStratification` stores exact
    % constant-stratification basis sets for $$N^2(z)=N_0^2$$. The class
    % implements the same `F`, `G`, and normalization contract as numerical
    % basis sets, without storing a solver reference.
    %
    % ```matlab
    % zDomain = [-5000 0];
    % N2 = @(z) (5.2e-3)^2*ones(size(z));
    % evp = IMInternalModes.waveModesAtWavenumber(N2=N2, zDomain=zDomain, k=1e-4);
    % basisSet = IMBasisSetConstantStratification(evp=evp, N0=5.2e-3, zDomain=[-5000 0]);
    % G = basisSet.G(linspace(-5000,0,128).');
    % ```
    %
    % - Topic: Create basis sets
    % - Topic: Evaluate basis sets
    % - Topic: Inspect basis sets
    % - Topic: Developer topics
    % - Declaration: classdef IMBasisSetConstantStratification < IMInternalModesBasis

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
        % Values are `"null"`, `"baroclinic"`, `"linear"`,
        % `"hyperbolic"`, or `"trig"`.
        %
        % - Topic: Inspect basis sets
        solutionTypes

        % True for a free-surface boundary branch.
        %
        % - Topic: Inspect basis sets
        isBoundaryMode
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
            % - Parameter options.evp: supported wave-mode or hydrostatic eigenvalue-problem descriptor
            % - Parameter options.N0: constant buoyancy frequency
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.nModes: number of retained modes
            % - Parameter options.normalization: active normalization rule; omitted uses the EVP default
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: exact constant-stratification basis set
            arguments
                options.evp = []
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.zDomain (1,2) double = [-1 0]
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.normalization = []
                options.metadata struct = struct()
            end

            zDomain = sort(options.zDomain);
            N2 = @(z) options.N0*options.N0*ones(size(z));
            evp = IMBasisSetConstantStratification.resolveEVP(options.evp, N2, zDomain);
            [f0, g] = IMBasisSetConstantStratification.physicalConstants(evp);
            [h, verticalWavenumbers, solutionTypes, isBoundaryMode, baroclinicNumbers, modeNumber] = ...
                IMBasisSetConstantStratification.solveSpectrum(evp, options.N0, zDomain, options.nModes, f0, g);
            eigenvalues = 1 ./ h;
            modeSelectionDiagnostics = struct();
            metadata = options.metadata;
            metadata.analyticalBasis = "constantStratification";

            self@IMInternalModesBasis(evp=evp, nativeModes=zeros(0,length(h)), ...
                eigenvalues=eigenvalues, h=h, modeNumber=modeNumber, modeSelectionDiagnostics=modeSelectionDiagnostics, ...
                normalization=options.normalization, metadata=metadata, zDomain=zDomain, ...
                N2=N2);
            self.N0 = options.N0;
            self.verticalWavenumbers = verticalWavenumbers;
            self.solutionTypes = solutionTypes;
            self.isBoundaryMode = isBoundaryMode;
            self.baroclinicNumbers = baroclinicNumbers;
        end

        function factors = normalizationFactors(self, normalization)
            % Return EVP-defined normalization factors.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: factors = normalizationFactors(basisSet,normalization)
            % - Parameter normalization: normalization rule name or enum value
            % - Returns factors: row vector of normalization factors
            arguments
                self IMBasisSetConstantStratification
                normalization = self.normalization
            end

            factors = normalizationFactors@IMBasisSet(self, normalization);
        end
    end

    methods (Hidden)
        function factor = maxAbsFactor(self, iMode, options)
            % Return an exact constant-stratification maximum-amplitude factor.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = maxAbsFactor(basisSet,iMode,options)
            % - Parameter iMode: mode index
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns factor: maximum absolute variable amplitude
            arguments
                self IMBasisSetConstantStratification
                iMode (1,1) double {mustBeInteger, mustBePositive}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            variable = string(options.variable);
            Lz = diff(self.zDomain);
            k_z = self.verticalWavenumbers(iMode);
            hMode = self.h(iMode);
            if ~self.isBoundaryMode(iMode)
                if self.evp.formulation == "F"
                    if variable == "G"
                        factor = abs(self.evp.g*k_z/(self.N0*self.N0));
                    else
                        factor = 1;
                    end
                elseif variable == "G"
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
                    if variable == "G"
                        factor = abs(sinh(k_z*Lz));
                    else
                        factor = abs(hMode*k_z*cosh(k_z*Lz));
                    end
                case "trig"
                    if variable == "G"
                        factor = abs(sin(k_z*Lz));
                    else
                        factor = abs(hMode*k_z);
                    end
                otherwise
                    factor = maxAbsFactor@IMInternalModesBasis(self, iMode, variable=variable);
            end
        end
    end

    methods (Access = protected)
        function values = rawVariable(self, variable, z)
            variable = string(variable);
            if variable ~= "G" && variable ~= "F"
                self.unsupported("evaluate " + variable);
            end

            z = z(:);
            s = z - self.zDomain(1);
            values = zeros(length(z), length(self.h));
            for iMode = 1:length(self.h)
                k_z = self.verticalWavenumbers(iMode);
                hMode = self.h(iMode);
                if self.isBoundaryMode(iMode)
                    [G, F] = self.rawBoundaryMode(self.solutionTypes(iMode), k_z, hMode, s);
                elseif self.evp.formulation == "F"
                    [G, F] = self.rawHydrostaticFMode(k_z, self.baroclinicNumbers(iMode), s);
                else
                    signValue = (-1)^self.baroclinicNumbers(iMode);
                    G = signValue*sin(k_z*s);
                    F = signValue*hMode*k_z*cos(k_z*s);
                end

                if variable == "G"
                    values(:,iMode) = G;
                else
                    values(:,iMode) = F;
                end
            end
        end
    end

    methods (Access = private)
        function [interiorIntegral, surfaceValue, bottomValue] = variableQuadraticPieces(self, variable, iMode)
            Lz = diff(self.zDomain);
            k_z = self.verticalWavenumbers(iMode);
            hMode = self.h(iMode);
            if ~self.isBoundaryMode(iMode)
                if self.evp.formulation == "F"
                    [interiorIntegral, surfaceValue, bottomValue] = self.hydrostaticFVariablePieces(variable, k_z, Lz);
                else
                    [interiorIntegral, surfaceValue, bottomValue] = ...
                        self.trigonometricVariablePieces(variable, k_z, hMode, Lz);
                end
                return;
            end

            switch string(self.solutionTypes(iMode))
                case "linear"
                    if variable == "G"
                        interiorIntegral = Lz*Lz*Lz/3;
                        surfaceValue = Lz;
                        bottomValue = 0;
                    else
                        interiorIntegral = Lz*Lz*Lz;
                        surfaceValue = Lz;
                        bottomValue = Lz;
                    end
                case "hyperbolic"
                    if variable == "G"
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
                        self.trigonometricVariablePieces(variable, k_z, hMode, Lz);
                otherwise
                    error("IMBasisSetConstantStratification:InvalidSolutionType", ...
                        "Unknown solution type ""%s"".", string(self.solutionTypes(iMode)));
            end
        end

        function [interiorIntegral, surfaceValue, bottomValue] = ...
                trigonometricVariablePieces(~, variable, k_z, hMode, Lz)
            if abs(k_z) <= eps
                interiorIntegral = 0;
                surfaceValue = 0;
                bottomValue = 0;
                return;
            end

            if variable == "G"
                interiorIntegral = Lz/2 - sin(2*k_z*Lz)/(4*k_z);
                surfaceValue = sin(k_z*Lz);
                bottomValue = 0;
            else
                interiorIntegral = hMode*hMode*k_z*k_z*(Lz/2 + sin(2*k_z*Lz)/(4*k_z));
                surfaceValue = hMode*k_z*cos(k_z*Lz);
                bottomValue = hMode*k_z;
            end
        end

        function [interiorIntegral, surfaceValue, bottomValue] = hydrostaticFVariablePieces(self, variable, k_z, Lz)
            if abs(k_z) <= eps
                if variable == "F"
                    interiorIntegral = Lz;
                    surfaceValue = 1;
                    bottomValue = 1;
                else
                    interiorIntegral = 0;
                    surfaceValue = 0;
                    bottomValue = 0;
                end
                return;
            end

            if variable == "F"
                interiorIntegral = Lz/2 + sin(2*k_z*Lz)/(4*k_z);
                surfaceValue = cos(k_z*Lz);
                bottomValue = 1;
            else
                scale = self.evp.g*k_z/(self.N0*self.N0);
                interiorIntegral = scale*scale*(Lz/2 - sin(2*k_z*Lz)/(4*k_z));
                surfaceValue = scale*sin(k_z*Lz);
                bottomValue = 0;
            end
        end

        function [G, F] = rawHydrostaticFMode(self, k_z, modeNumber, s)
            if abs(k_z) <= eps
                F = ones(size(s));
                G = zeros(size(s));
                return;
            end

            signValue = (-1)^modeNumber;
            F = signValue*cos(k_z*s);
            G = signValue*(self.evp.g*k_z/(self.N0*self.N0))*sin(k_z*s);
        end

        function [G, F] = rawBoundaryMode(self, solutionType, k_z, hMode, s)
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
        function evp = resolveEVP(evp, N2, zDomain)
            if isempty(evp)
                evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
                return;
            end
            if ~isa(evp, "IMInternalModes")
                error("IMBasisSetConstantStratification:InvalidEVP", ...
                    "The EVP must be an IMInternalModes instance.");
            end
            IMBasisSetConstantStratification.validateEVPDomain(evp, zDomain);
        end

        function validateEVPDomain(evp, zDomain)
            tolerance = 100*eps(max([1 abs(evp.zDomain) abs(zDomain)]));
            if max(abs(evp.zDomain - zDomain)) > tolerance
                error("IMBasisSetConstantStratification:DomainMismatch", ...
                    "The analytical basis zDomain must match evp.zDomain.");
            end
        end

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

        function [h, verticalWavenumbers, solutionTypes, isBoundaryMode, baroclinicNumbers, modeNumber] = ...
                solveSpectrum(evp, N0, zDomain, nModes, f0, g)
            [surfaceBoundary, bottomBoundary] = IMBasisSetConstantStratification.validateEVP(evp);
            D = diff(zDomain);
            evpName = string(evp.name);

            switch evpName
                case "waveModesAtWavenumber"
                    if ~isfield(evp.parameters, "k")
                        error("IMBasisSetConstantStratification:UnsupportedEVP", ...
                            "A fixed-wavenumber EVP must include parameters.k.");
                    end
                    k = evp.parameters.k;
                    [hBaroclinic, k_zBaroclinic] = IMBasisSetConstantStratification.baroclinicAtWavenumber( ...
                        k, N0, D, nModes, f0, g, surfaceBoundary);
                    if surfaceBoundary == "free"
                        [h0, k_z0, solutionType0] = IMBasisSetConstantStratification.surfaceBoundaryAtWavenumber( ...
                            k, N0, D, f0, g);
                        h = [h0 hBaroclinic(1:end-1)];
                        verticalWavenumbers = [k_z0 k_zBaroclinic(1:end-1)];
                        solutionTypes = IMBasisSetConstantStratification.freeSolutionTypes(solutionType0, nModes);
                        isBoundaryMode = [true false(1,nModes-1)];
                        baroclinicNumbers = [0 1:(nModes-1)];
                    else
                        h = hBaroclinic;
                        verticalWavenumbers = k_zBaroclinic;
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBoundaryMode = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "waveModesAtFrequency"
                    if ~isfield(evp.parameters, "omega")
                        error("IMBasisSetConstantStratification:UnsupportedEVP", ...
                            "A fixed-frequency EVP must include parameters.omega.");
                    end
                    omega = evp.parameters.omega;
                    if surfaceBoundary == "free"
                        [h0, k_z0, solutionType0] = IMBasisSetConstantStratification.surfaceBoundaryAtFrequency(omega, N0, D, g);
                        if omega < N0
                            [hBaroclinic, k_zBaroclinic] = IMBasisSetConstantStratification.freeSurfaceBaroclinicAtFrequency( ...
                                omega, N0, D, max(nModes - 1,0), g);
                        else
                            hBaroclinic = zeros(1,0);
                            k_zBaroclinic = zeros(1,0);
                        end
                        h = [h0 hBaroclinic];
                        verticalWavenumbers = [k_z0 k_zBaroclinic];
                        solutionTypes = IMBasisSetConstantStratification.freeSolutionTypes(solutionType0, length(h));
                        isBoundaryMode = [true false(1,length(h)-1)];
                        baroclinicNumbers = [0 1:(length(h)-1)];
                    else
                        if omega >= N0
                            error("IMBasisSetConstantStratification:UnsupportedFrequency", ...
                                "Rigid-surface fixed-frequency constant-stratification modes require omega < N0.");
                        end
                        [hBaroclinic, k_zBaroclinic] = IMBasisSetConstantStratification.baroclinicAtFrequency( ...
                            omega, N0, D, nModes, g);
                        h = hBaroclinic;
                        verticalWavenumbers = k_zBaroclinic;
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBoundaryMode = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "hydrostaticGModes"
                    if surfaceBoundary == "free"
                        [h0, k_z0, solutionType0] = IMBasisSetConstantStratification.surfaceBoundaryAtFrequency( ...
                            0, N0, D, g);
                        [hBaroclinic, k_zBaroclinic] = IMBasisSetConstantStratification.freeSurfaceBaroclinicAtFrequency( ...
                            0, N0, D, max(nModes - 1,0), g);
                        h = [h0 hBaroclinic];
                        verticalWavenumbers = [k_z0 k_zBaroclinic];
                        solutionTypes = IMBasisSetConstantStratification.freeSolutionTypes(solutionType0, length(h));
                        isBoundaryMode = [true false(1,length(h)-1)];
                        baroclinicNumbers = [0 1:(length(h)-1)];
                    else
                        [h, verticalWavenumbers] = IMBasisSetConstantStratification.baroclinicAtFrequency( ...
                            0, N0, D, nModes, g);
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBoundaryMode = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "hydrostaticFModes"
                    if surfaceBoundary ~= "rigid" || bottomBoundary ~= "rigid"
                        error("IMBasisSetConstantStratification:UnsupportedBoundary", ...
                            "Hydrostatic F constant-stratification modes currently require " + ...
                            "rigid surface and rigid bottom boundaries.");
                    end
                    interiorNumbers = 1:(nModes-1);
                    verticalWavenumbers = [0 interiorNumbers*pi/D];
                    h = [Inf (N0*N0)./(g*verticalWavenumbers(2:end).*verticalWavenumbers(2:end))];
                    solutionTypes = ["null" repmat("baroclinic",1,nModes-1)];
                    isBoundaryMode = false(1,nModes);
                    baroclinicNumbers = 0:(nModes-1);
                otherwise
                    error("IMBasisSetConstantStratification:UnsupportedEVP", ...
                        "Constant stratification supports fixed-wavenumber, fixed-frequency, " + ...
                        "hydrostatic G, and hydrostatic F EVPs.");
            end
            modeNumber = baroclinicNumbers;
            modeNumber(isBoundaryMode) = -1;
        end

        function [h, k_z] = baroclinicAtWavenumber(k, N0, D, nModes, f0, g, surfaceBoundary)
            k_z = (1:nModes)*pi/D;
            if surfaceBoundary == "free"
                for iMode = 1:nModes
                    f = @(xi) (xi + iMode*pi)*(N0*N0 - f0*f0)*D - g*(k*k*D*D + (xi + iMode*pi)*(xi + iMode*pi))*tan(xi);
                    k_z(iMode) = k_z(iMode) + fzero(f, 0)/D;
                end
            end
            h = (N0*N0 - f0*f0)./(g*(k*k + k_z.*k_z));
        end

        function solutionTypes = freeSolutionTypes(boundaryType, nModes)
            solutionTypes = repmat("baroclinic",1,nModes);
            solutionTypes(1) = boundaryType;
        end

        function [h, k_z] = baroclinicAtFrequency(omega, N0, D, nModes, g)
            k_z = (1:nModes)*pi/D;
            if omega >= N0
                error("IMBasisSetConstantStratification:UnsupportedFrequency", ...
                    "Interior fixed-frequency constant-stratification modes require omega < N0.");
            end
            h = (N0*N0 - omega*omega)./(g*k_z.*k_z);
        end

        function [h, k_z] = freeSurfaceBaroclinicAtFrequency(omega, N0, D, nModes, g)
            if omega >= N0
                error("IMBasisSetConstantStratification:UnsupportedFrequency", ...
                    "Free-surface interior fixed-frequency constant-stratification modes require omega < N0.");
            end
            k_z = (1:nModes)*pi/D;
            for iMode = 1:nModes
                f = @(xi) g*tan(xi)/(xi + iMode*pi) - (N0*N0 - omega*omega)*D/((xi + iMode*pi)^2);
                k_z(iMode) = k_z(iMode) + fzero(f, 0)/D;
            end
            h = (N0*N0 - omega*omega)./(g*k_z.*k_z);
        end

        function [h0, k_z, solutionType] = surfaceBoundaryAtWavenumber(k, N0, D, f0, g)
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

        function [h0, k_z, solutionType] = surfaceBoundaryAtFrequency(omega, N0, D, g)
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

        function [surfaceBoundary, bottomBoundary] = validateEVP(evp)
            surfaceBoundary = IMBasisSetConstantStratification.boundaryKind(evp.surfaceBoundary, evp.formulation);
            bottomBoundary = IMBasisSetConstantStratification.boundaryKind(evp.bottomBoundary, evp.formulation);
            if evp.formulation == "F"
                if evp.name ~= "hydrostaticFModes"
                    error("IMBasisSetConstantStratification:UnsupportedEVP", ...
                        "Constant stratification supports F-formulation EVPs only for hydrostatic F modes.");
                end
                if surfaceBoundary ~= "rigid" || bottomBoundary ~= "rigid"
                    error("IMBasisSetConstantStratification:UnsupportedBoundary", ...
                        "Hydrostatic F constant-stratification modes currently require " + ...
                        "rigid surface and rigid bottom boundaries.");
                end
                return;
            end
            if evp.formulation ~= "G"
                error("IMBasisSetConstantStratification:UnsupportedEVP", ...
                    "Constant stratification supports only F or G formulations.");
            end
            if bottomBoundary ~= "rigid"
                error("IMBasisSetConstantStratification:UnsupportedBoundary", ...
                    "Unsupported constant-stratification bottom boundary ""%s"".", bottomBoundary);
            end
            if ~ismember(surfaceBoundary, ["rigid", "free"])
                error("IMBasisSetConstantStratification:UnsupportedBoundary", ...
                    "Unsupported constant-stratification surface boundary ""%s"".", surfaceBoundary);
            end
        end

        function kind = boundaryKind(boundary, formulation)
            formulation = string(formulation);
            if formulation == "G"
                if IMBasisSetConstantStratification.isCondition(boundary, [1 0 0 0])
                    kind = "rigid";
                elseif IMBasisSetConstantStratification.isCondition(boundary, [0 1 1 0])
                    kind = "free";
                elseif IMBasisSetConstantStratification.isCondition(boundary, [0 1 0 0])
                    kind = "noSlip";
                else
                    kind = "custom";
                end
            else
                if IMBasisSetConstantStratification.isCondition(boundary, [0 1 0 0])
                    kind = "rigid";
                elseif IMBasisSetConstantStratification.isCondition(boundary, [1 0 0 0])
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
    end
end
