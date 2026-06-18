classdef IMConstantStratificationSolution < IMAnalyticalSolution
    % Analytical solution family for constant stratification.
    %
    % `IMConstantStratificationSolution` owns the closed-form formulas for
    % $$N^2(z)=N_0^2$$. It can create exact internal-mode bases for the
    % supported canonical internal-mode EVPs and exact surface or bottom SQG
    % boundary modes.
    %
    % ```matlab
    % solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
    % evp = IMInternalModes.hydrostaticGModes(N2=@(z) solution.N2(z), zDomain=solution.zDomain);
    % basisSet = solution.internalModes(evp,nModes=4);
    % ```
    %
    % - Topic: Create analytical solutions
    % - Topic: Compute internal modes
    % - Topic: Compute SQG modes
    % - Topic: Inspect analytical solutions
    % - Topic: Developer topics
    % - Declaration: classdef IMConstantStratificationSolution < IMAnalyticalSolution

    properties (SetAccess = private)
        % Constant buoyancy frequency $$N_0$$ in radians per second.
        %
        % - Topic: Inspect analytical solutions
        N0
    end

    methods
        function self = IMConstantStratificationSolution(options)
            % Create a constant-stratification analytical solution family.
            %
            % - Topic: Create analytical solutions
            % - Declaration: solution = IMConstantStratificationSolution(options)
            % - Parameter options.N0: constant buoyancy frequency
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Returns solution: analytical solution family
            arguments
                options.N0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 5.2e-3
                options.zDomain (1,2) double {mustBeReal, mustBeFinite} = [-1 0]
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
            end

            self@IMAnalyticalSolution(zDomain=options.zDomain, f0=options.f0, g=options.g);
            self.N0 = options.N0;
        end

        function values = N2(self, z)
            % Evaluate $$N^2(z)=N_0^2$$.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: values = N2(solution,z)
            % - Parameter z: physical coordinate
            % - Returns values: buoyancy frequency squared
            arguments
                self IMConstantStratificationSolution
                z double {mustBeReal, mustBeFinite}
            end

            values = self.N0*self.N0*ones(size(z));
        end

        function availability = internalModeAvailability(self, evp)
            % Report whether exact internal modes are available.
            %
            % - Topic: Compute internal modes
            % - Declaration: availability = internalModeAvailability(solution,evp)
            % - Parameter evp: internal-mode EVP
            % - Returns availability: availability report struct
            arguments
                self IMConstantStratificationSolution
                evp = []
            end

            try
                evp = self.resolveEVP(evp);
                IMConstantStratificationSolution.validateEVP(evp);
                availability = self.availabilityStruct(true, "internalModes", "constantStratification", "Constant-stratification formulas are available for this EVP.");
                availability.supportedNormalizations = IMConstantStratificationSolution.supportedNormalizations(evp);
            catch exception
                availability = self.availabilityStruct(false, "internalModes", "constantStratification", exception.message);
                availability.supportedNormalizations = string.empty(1,0);
            end
        end

        function basisSet = internalModes(self, evp, options)
            % Create an exact internal-mode basis.
            %
            % - Topic: Compute internal modes
            % - Declaration: basisSet = internalModes(solution,evp,options)
            % - Parameter evp: internal-mode EVP
            % - Parameter options.nModes: number of retained modes
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: exact analytical internal-mode basis
            arguments
                self IMConstantStratificationSolution
                evp = []
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.normalization = []
                options.metadata struct = struct()
            end

            evp = self.resolveEVP(evp);
            [f0, gValue] = IMConstantStratificationSolution.physicalConstants(evp);
            [h, verticalWavenumbers, solutionTypes, isBoundaryMode, baroclinicNumbers, modeNumber] = IMConstantStratificationSolution.solveSpectrum(evp, self.N0, self.zDomain, options.nModes, f0, gValue);
            modeData = struct();
            modeData.h = h;
            modeData.verticalWavenumbers = verticalWavenumbers;
            modeData.solutionTypes = solutionTypes;
            modeData.isBoundaryMode = isBoundaryMode;
            modeData.baroclinicNumbers = baroclinicNumbers;

            metadata = options.metadata;
            metadata.analyticalSolution = "constantStratification";
            metadata.verticalWavenumbers = verticalWavenumbers;
            metadata.solutionTypes = solutionTypes;
            metadata.isBoundaryMode = isBoundaryMode;
            metadata.baroclinicNumbers = baroclinicNumbers;

            basisSet = IMAnalyticalInternalModesBasis(solution=self, evp=evp, h=h, modeNumber=modeNumber, N2=@(z) self.N2(z), rawVariableFunction=@(basisSet,variable,z) self.rawVariable(modeData, evp, variable, z), rawUzFunction=@(basisSet,z) self.rawUz(modeData, evp, z), normalization=options.normalization, metadata=metadata);
        end

        function availability = sqgAvailability(self, options)
            % Report whether exact SQG modes are available.
            %
            % - Topic: Compute SQG modes
            % - Declaration: availability = sqgAvailability(solution,options)
            % - Parameter options.boundary: `"surface"` or `"bottom"`
            % - Returns availability: availability report struct
            arguments
                self IMConstantStratificationSolution
                options.boundary {mustBeTextScalar, mustBeMember(options.boundary, ["surface", "bottom"])} = "surface"
            end

            if abs(self.f0) <= eps(max(1,abs(self.f0)))
                availability = self.availabilityStruct(false, "sqgModes", "constantStratification", "SQG modes require nonzero f0.");
                return;
            end
            availability = self.availabilityStruct(true, "sqgModes", "constantStratification", "Constant-stratification SQG formulas are available.");
            availability.boundary = string(options.boundary);
            availability.supportedVariables = "psi";
            availability.supportedInnerProducts = string.empty(1,0);
            availability.supportedNormalizations = string.empty(1,0);
        end

        function basisSet = sqgModesAtWavenumber(self, k, options)
            % Create exact SQG boundary modes at fixed wavenumber.
            %
            % - Topic: Compute SQG modes
            % - Declaration: sqg = sqgModesAtWavenumber(solution,k,options)
            % - Parameter k: horizontal wavenumbers
            % - Parameter options.boundary: `"surface"` or `"bottom"`
            % - Parameter options.metadata: additional metadata
            % - Returns sqg: exact SQG basis
            arguments
                self IMConstantStratificationSolution
                k double {mustBeReal, mustBeFinite, mustBePositive}
                options.boundary {mustBeTextScalar, mustBeMember(options.boundary, ["surface", "bottom"])} = "surface"
                options.metadata struct = struct()
            end

            availability = self.sqgAvailability(boundary=options.boundary);
            if ~availability.isAvailable
                error("IMConstantStratificationSolution:UnavailableSQG", "%s", availability.reason);
            end
            metadata = options.metadata;
            metadata.analyticalSolution = "constantStratification";
            metadata.boundary = string(options.boundary);
            basisSet = IMAnalyticalSQGBasis(solution=self, k=k, boundary=options.boundary, N2=@(z) self.N2(z), psiFunction=@(z) self.sqgPsi(k, options.boundary, z), metadata=metadata);
        end

        function summarize(self, evp)
            % Print a readable solution-family summary.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: summarize(solution,evp)
            % - Parameter evp: optional internal-mode EVP
            arguments
                self IMConstantStratificationSolution
                evp = []
            end

            summarize@IMAnalyticalSolution(self, evp);
            fprintf("  N0: %g\n", self.N0);
        end
    end

    methods (Access = private)
        function evp = resolveEVP(self, evp)
            if isempty(evp)
                evp = IMInternalModes.hydrostaticGModes(N2=@(z) self.N2(z), zDomain=self.zDomain, f0=self.f0, g=self.g);
                return;
            end
            if ~isa(evp, "IMInternalModes")
                error("IMConstantStratificationSolution:InvalidEVP", "The EVP must be an IMInternalModes instance.");
            end
            IMConstantStratificationSolution.validateEVPDomain(evp, self.zDomain);
        end

        function values = rawVariable(self, modeData, evp, variable, z)
            variable = string(variable);
            if variable ~= "G" && variable ~= "F"
                error("IMConstantStratificationSolution:UnsupportedVariable", "Variable must be ""F"" or ""G"".");
            end

            z = z(:);
            s = z - self.zDomain(1);
            values = zeros(length(z), length(modeData.h));
            for iMode = 1:length(modeData.h)
                k_z = modeData.verticalWavenumbers(iMode);
                hMode = modeData.h(iMode);
                if modeData.isBoundaryMode(iMode)
                    [G, F] = self.rawBoundaryMode(modeData.solutionTypes(iMode), k_z, hMode, s);
                elseif evp.formulation == "F"
                    [G, F] = self.rawHydrostaticFMode(k_z, modeData.baroclinicNumbers(iMode), s, evp.g);
                else
                    signValue = (-1)^modeData.baroclinicNumbers(iMode);
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

        function values = rawUz(self, modeData, evp, z)
            z = z(:);
            s = z - self.zDomain(1);
            values = zeros(length(z), length(modeData.h));
            for iMode = 1:length(modeData.h)
                k_z = modeData.verticalWavenumbers(iMode);
                hMode = modeData.h(iMode);
                if modeData.isBoundaryMode(iMode)
                    [dGdz, dFdz] = self.rawBoundaryModeDerivatives(modeData.solutionTypes(iMode), k_z, hMode, s);
                elseif evp.formulation == "F"
                    [dGdz, dFdz] = self.rawHydrostaticFModeDerivatives(k_z, modeData.baroclinicNumbers(iMode), s, evp.g);
                else
                    signValue = (-1)^modeData.baroclinicNumbers(iMode);
                    dGdz = signValue*k_z*cos(k_z*s);
                    dFdz = -signValue*hMode*k_z*k_z*sin(k_z*s);
                end

                if evp.formulation == "G"
                    values(:,iMode) = dGdz;
                else
                    values(:,iMode) = dFdz;
                end
            end
        end

        function [G, F] = rawHydrostaticFMode(self, k_z, modeNumber, s, gValue)
            if abs(k_z) <= eps
                F = ones(size(s));
                G = zeros(size(s));
                return;
            end

            signValue = (-1)^modeNumber;
            F = signValue*cos(k_z*s);
            G = signValue*(gValue*k_z/(self.N0*self.N0))*sin(k_z*s);
        end

        function [dGdz, dFdz] = rawHydrostaticFModeDerivatives(self, k_z, modeNumber, s, gValue)
            if abs(k_z) <= eps
                dFdz = zeros(size(s));
                dGdz = zeros(size(s));
                return;
            end

            signValue = (-1)^modeNumber;
            dFdz = -signValue*k_z*sin(k_z*s);
            dGdz = signValue*(gValue*k_z*k_z/(self.N0*self.N0))*cos(k_z*s);
        end

        function [G, F] = rawBoundaryMode(self, solutionType, k_z, hMode, s)
            depth = diff(self.zDomain);
            switch string(solutionType)
                case "linear"
                    G = s;
                    F = depth*ones(size(s));
                case "hyperbolic"
                    G = sinh(k_z*s);
                    F = hMode*k_z*cosh(k_z*s);
                case "trig"
                    G = sin(k_z*s);
                    F = hMode*k_z*cos(k_z*s);
                otherwise
                    error("IMConstantStratificationSolution:InvalidSolutionType", "Unknown solution type ""%s"".", string(solutionType));
            end
        end

        function [dGdz, dFdz] = rawBoundaryModeDerivatives(~, solutionType, k_z, hMode, s)
            switch string(solutionType)
                case "linear"
                    dGdz = ones(size(s));
                    dFdz = zeros(size(s));
                case "hyperbolic"
                    dGdz = k_z*cosh(k_z*s);
                    dFdz = hMode*k_z*k_z*sinh(k_z*s);
                case "trig"
                    dGdz = k_z*cos(k_z*s);
                    dFdz = -hMode*k_z*k_z*sin(k_z*s);
                otherwise
                    error("IMConstantStratificationSolution:InvalidSolutionType", "Unknown solution type ""%s"".", string(solutionType));
            end
        end

        function values = sqgPsi(self, k, boundary, z)
            k = reshape(k,1,[]);
            z = z(:) - self.zDomain(2);
            depth = diff(self.zDomain);
            lambda = k*(self.N0/self.f0);
            denominator = k.*(1 - exp(-2*lambda*depth));
            if string(boundary) == "surface"
                numerator = exp(z.*lambda) + exp(-(z + 2*depth).*lambda);
                values = (1/self.N0)*numerator./denominator;
            else
                numerator = exp((z - depth).*lambda) + exp(-(z + depth).*lambda);
                values = -(1/self.N0)*numerator./denominator;
            end
        end
    end

    methods (Static, Access = private)
        function validateEVPDomain(evp, zDomain)
            tolerance = 100*eps(max([1 abs(evp.zDomain) abs(zDomain)]));
            if max(abs(evp.zDomain - zDomain)) > tolerance
                error("IMConstantStratificationSolution:DomainMismatch", "The analytical solution zDomain must match evp.zDomain.");
            end
        end

        function [f0, gValue] = physicalConstants(evp)
            f0 = evp.f0;
            gValue = evp.g;
            if ~(isscalar(f0) && isfinite(f0))
                error("IMConstantStratificationSolution:InvalidCoriolis", "The Coriolis parameter must be finite.");
            end
            if ~(isscalar(gValue) && isfinite(gValue) && gValue > 0)
                error("IMConstantStratificationSolution:InvalidGravity", "The gravitational acceleration must be positive.");
            end
        end

        function normalizations = supportedNormalizations(evp)
            normalizations = ["unity", "uMax", "wMax", "surfacePressure"];
            if evp.modeFamily == "hydrostatic"
                normalizations(end+1) = "geostrophic";
            end
            if string(evp.name) == "waveModesAtWavenumber"
                normalizations(end+1) = "kConstant";
            end
        end

        function [h, verticalWavenumbers, solutionTypes, isBoundaryMode, baroclinicNumbers, modeNumber] = solveSpectrum(evp, N0, zDomain, nModes, f0, gValue)
            [surfaceBoundary, bottomBoundary] = IMConstantStratificationSolution.validateEVP(evp);
            depth = diff(zDomain);
            evpName = string(evp.name);

            switch evpName
                case "waveModesAtWavenumber"
                    if ~isfield(evp.parameters, "k")
                        error("IMConstantStratificationSolution:UnsupportedEVP", "A fixed-wavenumber EVP must include parameters.k.");
                    end
                    k = evp.parameters.k;
                    [hBaroclinic, k_zBaroclinic] = IMConstantStratificationSolution.baroclinicAtWavenumber(k, N0, depth, nModes, f0, gValue, surfaceBoundary);
                    if surfaceBoundary == "free"
                        [h0, k_z0, solutionType0] = IMConstantStratificationSolution.surfaceBoundaryAtWavenumber(k, N0, depth, f0, gValue);
                        h = [h0 hBaroclinic(1:end-1)];
                        verticalWavenumbers = [k_z0 k_zBaroclinic(1:end-1)];
                        solutionTypes = IMConstantStratificationSolution.freeSolutionTypes(solutionType0, nModes);
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
                        error("IMConstantStratificationSolution:UnsupportedEVP", "A fixed-frequency EVP must include parameters.omega.");
                    end
                    omega = evp.parameters.omega;
                    if surfaceBoundary == "free"
                        [h0, k_z0, solutionType0] = IMConstantStratificationSolution.surfaceBoundaryAtFrequency(omega, N0, depth, gValue);
                        if omega < N0
                            [hBaroclinic, k_zBaroclinic] = IMConstantStratificationSolution.freeSurfaceBaroclinicAtFrequency(omega, N0, depth, max(nModes - 1,0), gValue);
                        else
                            hBaroclinic = zeros(1,0);
                            k_zBaroclinic = zeros(1,0);
                        end
                        h = [h0 hBaroclinic];
                        verticalWavenumbers = [k_z0 k_zBaroclinic];
                        solutionTypes = IMConstantStratificationSolution.freeSolutionTypes(solutionType0, length(h));
                        isBoundaryMode = [true false(1,length(h)-1)];
                        baroclinicNumbers = [0 1:(length(h)-1)];
                    else
                        if omega >= N0
                            error("IMConstantStratificationSolution:UnsupportedFrequency", "Rigid-surface fixed-frequency constant-stratification modes require omega < N0.");
                        end
                        [hBaroclinic, k_zBaroclinic] = IMConstantStratificationSolution.baroclinicAtFrequency(omega, N0, depth, nModes, gValue);
                        h = hBaroclinic;
                        verticalWavenumbers = k_zBaroclinic;
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBoundaryMode = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "hydrostaticGModes"
                    if surfaceBoundary == "free"
                        [h0, k_z0, solutionType0] = IMConstantStratificationSolution.surfaceBoundaryAtFrequency(0, N0, depth, gValue);
                        [hBaroclinic, k_zBaroclinic] = IMConstantStratificationSolution.freeSurfaceBaroclinicAtFrequency(0, N0, depth, max(nModes - 1,0), gValue);
                        h = [h0 hBaroclinic];
                        verticalWavenumbers = [k_z0 k_zBaroclinic];
                        solutionTypes = IMConstantStratificationSolution.freeSolutionTypes(solutionType0, length(h));
                        isBoundaryMode = [true false(1,length(h)-1)];
                        baroclinicNumbers = [0 1:(length(h)-1)];
                    else
                        [h, verticalWavenumbers] = IMConstantStratificationSolution.baroclinicAtFrequency(0, N0, depth, nModes, gValue);
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBoundaryMode = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "hydrostaticFModes"
                    if surfaceBoundary ~= "rigid" || bottomBoundary ~= "rigid"
                        error("IMConstantStratificationSolution:UnsupportedBoundary", "Hydrostatic F constant-stratification modes currently require rigid surface and rigid bottom boundaries.");
                    end
                    interiorNumbers = 1:(nModes-1);
                    verticalWavenumbers = [0 interiorNumbers*pi/depth];
                    h = [Inf (N0*N0)./(gValue*verticalWavenumbers(2:end).*verticalWavenumbers(2:end))];
                    solutionTypes = ["null" repmat("baroclinic",1,nModes-1)];
                    isBoundaryMode = false(1,nModes);
                    baroclinicNumbers = 0:(nModes-1);
                otherwise
                    error("IMConstantStratificationSolution:UnsupportedEVP", "Constant stratification supports fixed-wavenumber, fixed-frequency, hydrostatic G, and hydrostatic F EVPs.");
            end
            modeNumber = baroclinicNumbers;
            modeNumber(isBoundaryMode) = -1;
        end

        function [h, k_z] = baroclinicAtWavenumber(k, N0, depth, nModes, f0, gValue, surfaceBoundary)
            k_z = (1:nModes)*pi/depth;
            if surfaceBoundary == "free"
                for iMode = 1:nModes
                    residual = @(xi) (xi + iMode*pi)*(N0*N0 - f0*f0)*depth - gValue*(k*k*depth*depth + (xi + iMode*pi)*(xi + iMode*pi))*tan(xi);
                    k_z(iMode) = k_z(iMode) + fzero(residual, 0)/depth;
                end
            end
            h = (N0*N0 - f0*f0)./(gValue*(k*k + k_z.*k_z));
        end

        function solutionTypes = freeSolutionTypes(boundaryType, nModes)
            solutionTypes = repmat("baroclinic",1,nModes);
            solutionTypes(1) = boundaryType;
        end

        function [h, k_z] = baroclinicAtFrequency(omega, N0, depth, nModes, gValue)
            k_z = (1:nModes)*pi/depth;
            if omega >= N0
                error("IMConstantStratificationSolution:UnsupportedFrequency", "Interior fixed-frequency constant-stratification modes require omega < N0.");
            end
            h = (N0*N0 - omega*omega)./(gValue*k_z.*k_z);
        end

        function [h, k_z] = freeSurfaceBaroclinicAtFrequency(omega, N0, depth, nModes, gValue)
            if omega >= N0
                error("IMConstantStratificationSolution:UnsupportedFrequency", "Free-surface interior fixed-frequency constant-stratification modes require omega < N0.");
            end
            k_z = (1:nModes)*pi/depth;
            for iMode = 1:nModes
                residual = @(xi) gValue*tan(xi)/(xi + iMode*pi) - (N0*N0 - omega*omega)*depth/((xi + iMode*pi)^2);
                k_z(iMode) = k_z(iMode) + fzero(residual, 0)/depth;
            end
            h = (N0*N0 - omega*omega)./(gValue*k_z.*k_z);
        end

        function [h0, k_z, solutionType] = surfaceBoundaryAtWavenumber(k, N0, depth, f0, gValue)
            kStar = sqrt((N0*N0 - f0*f0)/(gValue*depth));
            if abs(k - kStar)/kStar < 1e-6
                solutionType = "linear";
                h0 = depth;
                k_z = 0;
            elseif k > kStar
                solutionType = "hyperbolic";
                residual = @(q) depth*(N0*N0 - f0*f0) - (1./q).*(gValue*(k*k*depth*depth - q.*q)).*tanh(q);
                kInitial = sqrt(k*k*depth*depth - depth*(N0*N0 - f0*f0)/gValue);
                k_z = fzero(residual, kInitial)/depth;
                h0 = (N0*N0 - f0*f0)/(gValue*(k*k - k_z*k_z));
            else
                solutionType = "trig";
                residual = @(q) depth*(N0*N0 - f0*f0) - (1./q).*(gValue*(k*k*depth*depth + q.*q)).*tan(q);
                kInitial = sqrt(-k*k*depth*depth + depth*(N0*N0 - f0*f0)/gValue);
                k_z = fzero(residual, kInitial)/depth;
                h0 = (N0*N0 - f0*f0)/(gValue*(k*k + k_z*k_z));
            end
        end

        function [h0, k_z, solutionType] = surfaceBoundaryAtFrequency(omega, N0, depth, gValue)
            if abs(omega - N0)/N0 < 1e-6
                solutionType = "linear";
                h0 = depth;
                k_z = 0;
            elseif omega > N0
                solutionType = "hyperbolic";
                residual = @(q) depth*(omega*omega - N0*N0) - gValue*q.*tanh(q);
                k_z = fzero(residual, sqrt(depth*(omega*omega - N0*N0)/gValue))/depth;
                h0 = (omega*omega - N0*N0)/(gValue*k_z*k_z);
            else
                solutionType = "trig";
                residual = @(q) depth*(N0*N0 - omega*omega) - gValue*q.*tan(q);
                k_z = fzero(residual, sqrt(depth*(N0*N0 - omega*omega)/gValue))/depth;
                h0 = (N0*N0 - omega*omega)/(gValue*k_z*k_z);
            end
        end

        function [surfaceBoundary, bottomBoundary] = validateEVP(evp)
            surfaceBoundary = IMConstantStratificationSolution.boundaryKind(evp.surfaceBoundary, evp.formulation);
            bottomBoundary = IMConstantStratificationSolution.boundaryKind(evp.bottomBoundary, evp.formulation);
            if evp.formulation == "F"
                if evp.name ~= "hydrostaticFModes"
                    error("IMConstantStratificationSolution:UnsupportedEVP", "Constant stratification supports F-formulation EVPs only for hydrostatic F modes.");
                end
                if surfaceBoundary ~= "rigid" || bottomBoundary ~= "rigid"
                    error("IMConstantStratificationSolution:UnsupportedBoundary", "Hydrostatic F constant-stratification modes currently require rigid surface and rigid bottom boundaries.");
                end
                return;
            end
            if evp.formulation ~= "G"
                error("IMConstantStratificationSolution:UnsupportedEVP", "Constant stratification supports only F or G formulations.");
            end
            if bottomBoundary ~= "rigid"
                error("IMConstantStratificationSolution:UnsupportedBoundary", "Unsupported constant-stratification bottom boundary ""%s"".", bottomBoundary);
            end
            if ~ismember(surfaceBoundary, ["rigid", "free"])
                error("IMConstantStratificationSolution:UnsupportedBoundary", "Unsupported constant-stratification surface boundary ""%s"".", surfaceBoundary);
            end
        end

        function kind = boundaryKind(boundary, formulation)
            formulation = string(formulation);
            if formulation == "G"
                if IMConstantStratificationSolution.isCondition(boundary, [1 0 0 0])
                    kind = "rigid";
                elseif IMConstantStratificationSolution.isCondition(boundary, [0 1 1 0])
                    kind = "free";
                elseif IMConstantStratificationSolution.isCondition(boundary, [0 1 0 0])
                    kind = "noSlip";
                else
                    kind = "custom";
                end
            else
                if IMConstantStratificationSolution.isCondition(boundary, [0 1 0 0])
                    kind = "rigid";
                elseif IMConstantStratificationSolution.isCondition(boundary, [1 0 0 0])
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
