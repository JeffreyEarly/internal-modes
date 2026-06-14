classdef IMBasisSet
    % Store solved scalar canonical EVP modes.
    %
    % `IMBasisSet` stores native coefficient columns returned by a solver
    % and evaluates the solved scalar variable `u` and its derivative. Modal
    % normalization is applied lazily when values or Gram matrices are
    % requested. For each mode,
    % $$u_j^{\mathrm{out}}(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
    % where the scale factor $$s_j$$ is supplied by the active
    % normalization rule.
    %
    % ```matlab
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet.normalization = Normalization.unity;
    % u = basisSet.u(z);
    % factors = basisSet.normalizationFactors(Normalization.unity);
    % ```
    %
    % - Topic: Create basis sets
    % - Topic: Evaluate basis sets
    % - Topic: Analyze Gram matrices
    % - Topic: Inspect basis sets
    % - Topic: Developer topics
    % - Declaration: classdef IMBasisSet

    properties
        % Active modal normalization.
        %
        % This value selects the normalization rule used by `u`, `uz`,
        % `evaluate`, and Gram-matrix methods. Passing `normalization=...`
        % to an evaluation method overrides this property for that call.
        %
        % - Topic: Evaluate basis sets
        normalization = Normalization.unity
    end

    properties (SetAccess = protected)
        % Solver that created the native modes.
        %
        % - Topic: Inspect basis sets
        solver

        % EVP descriptor that was solved.
        %
        % - Topic: Inspect basis sets
        evp

        % Native mode columns.
        %
        % - Topic: Inspect basis sets
        nativeModes

        % Retained eigenvalues.
        %
        % The retained values are the discrete eigenvalues
        % $$\lambda_j$$ from the assembled canonical problem.
        %
        % - Topic: Inspect basis sets
        eigenvalues

        % Physical mode numbers.
        %
        % - Topic: Inspect basis sets
        modeNumber

        % Diagnostic index information from mode selection.
        %
        % - Topic: Inspect basis sets
        index

        % Additional metadata.
        %
        % `metadata` stores creation or diagnostic information associated
        % with this solved basis set. It is not interpreted by `IMBasisSet`;
        % EVP coefficient parameters remain available through
        % `basisSet.evp.parameters`.
        %
        % - Topic: Inspect basis sets
        metadata

        % Physical vertical domain.
        %
        % - Topic: Inspect basis sets
        zDomain
    end

    properties (Access = protected)
        % Buoyancy frequency squared function.
        %
        % - Topic: Developer topics
        % - Developer: true
        N2Function
    end

    methods
        function self = IMBasisSet(options)
            % Create a solved scalar basis set.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet(options)
            % - Parameter options.solver: solver reference
            % - Parameter options.evp: canonical EVP descriptor
            % - Parameter options.nativeModes: native mode columns
            % - Parameter options.eigenvalues: retained eigenvalues
            % - Parameter options.modeNumber: physical mode numbers
            % - Parameter options.index: selection diagnostics
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.N2Function: buoyancy frequency squared function
            % - Returns basisSet: solved scalar basis set
            arguments
                options.solver = []
                options.evp = []
                options.nativeModes double = zeros(0,0)
                options.eigenvalues double = zeros(1,0)
                options.modeNumber double = []
                options.index struct = struct()
                options.normalization = []
                options.metadata struct = struct()
                options.zDomain (1,2) double = [NaN NaN]
                options.N2Function = []
            end

            self.solver = options.solver;
            self.evp = options.evp;
            self.nativeModes = options.nativeModes;
            self.eigenvalues = reshape(options.eigenvalues,1,[]);
            nModes = max(size(options.nativeModes,2), length(self.eigenvalues));
            self.modeNumber = IMBasisSet.resolveModeNumber(options.modeNumber, nModes);
            self.index = options.index;
            self.normalization = IMBasisSet.resolveDefaultNormalization(options.normalization, options.evp);
            self.metadata = options.metadata;
            self.zDomain = options.zDomain;
            self.N2Function = options.N2Function;

            if any(isnan(self.zDomain)) && ~isempty(self.evp)
                self.zDomain = self.evp.zDomain;
            elseif ~isempty(self.solver)
                if any(isnan(self.zDomain))
                    self.zDomain = self.solver.zDomain;
                end
            end
            if isempty(self.N2Function) && isa(self.evp, "IMInternalModes")
                self.N2Function = self.evp.N2;
            end
        end

        function values = u(self, z, options)
            % Evaluate solved scalar modes.
            %
            % Returned columns are normalized as
            % $$u_j(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
            % where $$s_j$$ comes from `normalizationFactors`.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = u(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization to apply
            % - Returns values: scalar mode values
            arguments
                self IMBasisSet
                z (:,1) double
                options.normalization = self.normalization
            end

            values = self.rawU(z) ./ self.normalizationFactors(options.normalization);
        end

        function values = uz(self, z, options)
            % Evaluate solved scalar vertical derivatives.
            %
            % Derivatives are scaled by the same modal factors used for
            % `u`, so $$u'_j(z)=u_j^{\mathrm{raw}\prime}(z)/s_j$$.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = uz(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization to apply
            % - Returns values: derivative mode values
            arguments
                self IMBasisSet
                z (:,1) double
                options.normalization = self.normalization
            end

            values = self.rawUz(z) ./ self.normalizationFactors(options.normalization);
        end

        function values = evaluate(self, variable, z, options)
            % Evaluate the scalar variable.
            %
            % `evaluate("u",z)` is equivalent to `u(z)` and accepts the
            % same `normalization` override.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = evaluate(basisSet,variable,z,options)
            % - Parameter variable: `"u"`
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization to apply
            % - Returns values: scalar mode values
            arguments
                self IMBasisSet
                variable {mustBeTextScalar}
                z (:,1) double
                options.normalization = self.normalization
            end

            if string(variable) ~= "u"
                error("IMBasisSet:InvalidVariable", ...
                    "Canonical scalar basis sets evaluate variable ""u"".");
            end
            values = self.u(z, normalization=options.normalization);
        end

        function values = evaluateAll(self, z)
            % Evaluate all scalar fields.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = evaluateAll(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: struct containing `u` and `uz`
            values.u = self.u(z);
            values.uz = self.uz(z);
        end

        function factors = normalizationFactors(self, normalization)
            % Return factors for a normalization convention.
            %
            % For a requested convention this returns the row vector
            % $$s_j$$ used by evaluation methods:
            % $$u_j^{\mathrm{out}}=u_j^{\mathrm{raw}}/s_j.$$
            % The default scalar `unity` convention uses
            % $$s_j=\sqrt{|\langle u_j,u_j\rangle|}$$ with the EVP inner
            % product.
            %
            % ```matlab
            % factors = basisSet.normalizationFactors(Normalization.unity);
            % uUnity = basisSet.u(z,normalization=Normalization.unity);
            % ```
            %
            % - Topic: Evaluate basis sets
            % - Declaration: factors = normalizationFactors(basisSet,normalization)
            % - Parameter normalization: normalization convention
            % - Returns factors: row vector of positive scale factors
            name = char(self.normalizationName(normalization));
            if isempty(self.evp) || ~isfield(self.evp.normalizations, name)
                error("IMBasisSet:UnsupportedNormalization", ...
                    "The EVP does not define a ""%s"" normalization.", name);
            end
            normalizeMode = self.evp.normalizations.(name);
            nModes = max(length(self.eigenvalues), size(self.nativeModes,2));
            factors = zeros(1,nModes);
            for iMode = 1:nModes
                factors(iMode) = normalizeMode(self, iMode);
            end
            factors = abs(factors);
            factors(factors == 0 | ~isfinite(factors)) = 1;
        end

        function modes = normalizedNativeModes(self, normalization)
            % Return native modes scaled by a normalization.
            %
            % The native coefficient columns are divided by the same factors
            % returned by `normalizationFactors`.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: modes = normalizedNativeModes(basisSet,normalization)
            % - Parameter normalization: normalization convention
            % - Returns modes: scaled native mode columns
            factors = self.normalizationFactors(normalization);
            modes = self.nativeModes ./ factors;
        end

        function gram = gramMatrix(self, variable)
            % Return the full-domain scalar Gram matrix.
            %
            % For normalized scalar modes, entries are
            % $$M_{ij}=\int r u_i u_j\,dz+
            % \sum_\ell \gamma_\ell L_\ell[u_i]L_\ell[u_j].$$
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = gramMatrix(basisSet,variable)
            % - Parameter variable: optional scalar variable name, `"u"`
            % - Returns gram: scalar Gram matrix
            arguments
                self IMBasisSet
                variable {mustBeTextScalar} = "u"
            end

            gram = self.partialGramMatrix(variable, self.zDomain(1), self.zDomain(2));
        end

        function gram = partialGramMatrix(self, variable, zMin, zMax)
            % Return a partial-domain scalar Gram matrix.
            %
            % Interior integrals are restricted to `[zMin,zMax]`. Endpoint
            % metric terms are included only when the requested interval
            % contains the corresponding physical endpoint.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = partialGramMatrix(basisSet,variable,zMin,zMax)
            % - Parameter variable: scalar variable name, `"u"`
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns gram: scalar Gram matrix
            if string(variable) ~= "u"
                error("IMBasisSet:InvalidVariable", ...
                    "Canonical scalar basis sets use variable ""u"".");
            end
            gram = self.scalarGramMatrix([zMin zMax], true);
        end

        function windowModes = partialWindowModes(self, variable, zMin, zMax)
            % Diagonalize a partial scalar Gram matrix.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: windowModes = partialWindowModes(basisSet,variable,zMin,zMax)
            % - Parameter variable: scalar variable name
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns windowModes: window-mode decomposition
            gram = self.partialGramMatrix(variable, zMin, zMax);
            gram = 0.5*(gram + gram.');
            [R, D] = eig(gram);
            [eigenvalues, sortIndex] = sort(diag(D), "descend");
            windowModes.rotation = R(:,sortIndex);
            windowModes.eigenvalues = eigenvalues(:).';
            windowModes.gramMatrix = gram;
        end

        function spectrum = spectrum(self, coefficients, options)
            % Compute a scalar modal spectrum.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = spectrum(basisSet,coefficients,options)
            % - Parameter coefficients: modal coefficients
            % - Parameter options.variable: scalar variable name
            % - Returns spectrum: modal spectrum
            arguments
                self IMBasisSet
                coefficients (:,1) double
                options.variable {mustBeTextScalar} = "u"
            end

            spectrum = self.crossSpectrum(coefficients, coefficients, variable=options.variable);
        end

        function spectrum = crossSpectrum(self, coefficientsA, coefficientsB, options)
            % Compute a scalar modal cross-spectrum.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB,options)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Parameter options.variable: scalar variable name
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMBasisSet
                coefficientsA (:,1) double
                coefficientsB (:,1) double
                options.variable {mustBeTextScalar} = "u"
            end

            gram = self.gramMatrix(options.variable);
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
        end

        function values = N2(self, z)
            % Evaluate buoyancy frequency squared.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = N2(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: buoyancy frequency squared
            if isempty(self.N2Function)
                error("IMBasisSet:UnsupportedOperation", ...
                    "The basis set does not have an N2 function.");
            end
            values = self.N2Function(z);
        end

        function values = dzLogN2(self, z)
            % Evaluate the vertical derivative of `log(N2)`.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = dzLogN2(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: derivative values
            if isempty(self.N2Function)
                error("IMBasisSet:UnsupportedOperation", ...
                    "The basis set does not have an N2 function.");
            end

            z = z(:);
            if length(z) == 1
                scale = max(1,abs(z));
                dz = sqrt(eps)*scale;
                values = (log(self.N2Function(z + dz)) - log(self.N2Function(z - dz)))/(2*dz);
            else
                values = gradient(log(self.N2Function(z)), z);
            end
        end
    end

    methods (Hidden)
        function self = orientModeSigns(self)
            % Orient scalar modes with a deterministic sign convention.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = orientModeSigns(basisSet)
            % - Developer: true
            if isempty(self.nativeModes)
                return;
            end
            z = self.innerProductGrid(self.zDomain);
            values = self.rawU(z);
            signs = ones(1,size(values,2));
            for iMode = 1:size(values,2)
                [~, index] = max(abs(values(:,iMode)));
                if values(index,iMode) < 0
                    signs(iMode) = -1;
                end
            end
            self.nativeModes = self.nativeModes .* signs;
        end

        function factor = innerProductNormFactor(self, varargin)
            % Return the scalar inner-product norm factor.
            %
            % This is the raw factor
            % $$s_j=\sqrt{|\langle u_j,u_j\rangle|}$$ computed before any
            % modal normalization has been applied.
            %
            % - Topic: Developer topics
            % - Declaration: factor = innerProductNormFactor(basisSet,iMode)
            % - Developer: true
            if nargin == 2
                iMode = varargin{1};
            elseif nargin == 3 && string(varargin{1}) == "u"
                iMode = varargin{2};
            else
                error("IMBasisSet:InvalidVariable", ...
                    "Canonical scalar norm factors use variable ""u"".");
            end
            gram = self.scalarGramMatrix(self.zDomain, false);
            factor = sqrt(abs(gram(iMode,iMode)));
        end

        function factor = maxAbsFactor(self, varargin)
            % Return the maximum scalar amplitude.
            %
            % This is $$s_j=\max_z |u_j^{\mathrm{raw}}(z)|$$ on the
            % basis-set integration grid.
            %
            % - Topic: Developer topics
            % - Declaration: factor = maxAbsFactor(basisSet,iMode)
            % - Developer: true
            if nargin == 2
                iMode = varargin{1};
            elseif nargin == 3 && string(varargin{1}) == "u"
                iMode = varargin{2};
            else
                error("IMBasisSet:InvalidVariable", ...
                    "Canonical scalar maximum factors use variable ""u"".");
            end
            z = self.integrationGrid(self.zDomain);
            values = self.rawU(z);
            factor = max(abs(values(:,iMode)));
        end
    end

    methods (Access = protected)
        function values = rawU(self, z)
            self.requireSolver("evaluate solved scalar modes");
            values = self.solver.evaluateNativeModes(self.nativeModes, z);
        end

        function values = rawUz(self, z)
            self.requireSolver("evaluate scalar derivatives");
            values = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
        end

        function z = integrationGrid(self, zBounds)
            if ~isempty(self.solver)
                nGrid = max(256, 4*self.solver.nEVP);
            else
                nGrid = max(256, 4*max(1,length(self.eigenvalues)));
            end
            z = linspace(min(zBounds), max(zBounds), nGrid).';
        end

        function z = innerProductGrid(self, zBounds)
            if ~isempty(self.solver)
                z = self.solver.innerProductGrid(zBounds);
            else
                z = self.integrationGrid(zBounds);
            end
        end

        function value = integrateInnerProduct(self, z, integrand, zBounds)
            if ~isempty(self.solver)
                value = self.solver.integrateInnerProduct(z, integrand, zBounds);
            else
                value = trapz(z, integrand);
            end
        end

        function gram = scalarGramMatrix(self, zBounds, useNormalized)
            z = self.innerProductGrid(zBounds);
            if useNormalized
                values = self.u(z);
            else
                values = self.rawU(z);
            end
            spec = self.evp.innerProduct("u");
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, self.context());
            if isscalar(weight)
                weight = weight*ones(size(z));
            end
            gram = zeros(size(values,2), size(values,2));
            for iMode = 1:size(values,2)
                for jMode = iMode:size(values,2)
                    integrand = weight(:).*values(:,iMode).*values(:,jMode);
                    value = self.integrateInnerProduct(z, integrand, zBounds);
                    value = value + self.endpointMetricValue([spec.surfaceWeights; spec.bottomWeights], ...
                        iMode, jMode, useNormalized, zBounds);
                    gram(iMode,jMode) = value;
                    gram(jMode,iMode) = value;
                end
            end
        end

        function value = endpointMetricValue(self, weights, iMode, jMode, useNormalized, zBounds)
            value = 0;
            for iWeight = 1:numel(weights)
                weight = weights(iWeight);
                zEndpoint = self.endpointZ(weight.location);
                if ~self.boundsIncludeEndpoint(zBounds, zEndpoint)
                    continue;
                end
                left = self.endpointMetricFactor(weight, zEndpoint, useNormalized);
                value = value + weight.coefficient*left(iMode)*left(jMode);
            end
        end

        function values = endpointMetricFactor(self, weight, zEndpoint, useNormalized)
            if useNormalized
                uValues = self.u(zEndpoint);
                uzValues = self.uz(zEndpoint);
            else
                uValues = self.rawU(zEndpoint);
                uzValues = self.rawUz(zEndpoint);
            end
            pValue = IMEigenvalueProblem.evaluateCoefficient(self.evp.p, zEndpoint, self.context());
            values = weight.c*uValues - weight.d*pValue*uzValues;
        end

        function zEndpoint = endpointZ(self, location)
            switch string(location)
                case "surface"
                    zEndpoint = self.zDomain(2);
                case "bottom"
                    zEndpoint = self.zDomain(1);
                otherwise
                    error("IMBasisSet:InvalidBoundaryLocation", ...
                        "Boundary location must be ""surface"" or ""bottom"".");
            end
        end

        function context = context(self)
            context.N2 = @(z) self.N2(z);
            context.dzLogN2 = @(z) self.dzLogN2(z);
            context.zDomain = self.zDomain;
            context.coordinateKind = "basisSet";
            if ~isempty(self.solver)
                solverContext = self.solver.context();
                contextFields = fieldnames(solverContext);
                for iField = 1:numel(contextFields)
                    context.(contextFields{iField}) = solverContext.(contextFields{iField});
                end
            end
            if ~isempty(self.evp)
                parameterFields = fieldnames(self.evp.parameters);
                for iField = 1:numel(parameterFields)
                    fieldName = parameterFields{iField};
                    context.(fieldName) = self.evp.parameters.(fieldName);
                end
            end
        end

        function isIncluded = boundsIncludeEndpoint(~, zBounds, endpoint)
            tolerance = 100*eps(max(1,max(abs([zBounds(:); endpoint]))));
            isIncluded = abs(min(zBounds) - endpoint) <= tolerance || abs(max(zBounds) - endpoint) <= tolerance;
        end

        function name = normalizationName(~, normalization)
            name = erase(string(normalization), "Normalization.");
        end

        function requireSolver(self, operationName)
            if isempty(self.solver)
                error("IMBasisSet:UnsupportedOperation", ...
                    "The basis set does not have a solver reference and cannot perform %s.", operationName);
            end
        end

        function unsupported(~, operationName)
            error("IMBasisSet:UnsupportedOperation", ...
                "IMBasisSet does not support %s for this basis.", operationName);
        end
    end

    methods (Static)
        function basisSet = constantStratification(options)
            % Create an analytical constant-stratification basis set.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet.constantStratification(options)
            % - Parameter options.evp: internal-mode EVP descriptor
            % - Parameter options.N0: constant buoyancy frequency
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.nModes: number of modes
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: analytical constant-stratification basis set
            arguments
                options.evp = []
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.zDomain (1,2) double = [-1 0]
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.normalization = []
                options.metadata struct = struct()
            end

            basisSet = IMBasisSetConstantStratification(evp=options.evp, N0=options.N0, ...
                zDomain=options.zDomain, nModes=options.nModes, normalization=options.normalization, ...
                metadata=options.metadata);
        end

        function basisSet = exponentialStratification(options)
            % Create an analytical exponential-stratification basis set.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet.exponentialStratification(options)
            % - Parameter options.evp: internal-mode EVP descriptor
            % - Parameter options.N0: surface buoyancy frequency
            % - Parameter options.b: exponential e-folding depth
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.nModes: number of modes
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: analytical exponential-stratification basis set
            arguments
                options.evp = []
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.b (1,1) double {mustBePositive} = 1300
                options.zDomain (1,2) double = [-1 0]
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.normalization = []
                options.metadata struct = struct()
            end

            basisSet = IMBasisSetExponentialStratification(evp=options.evp, N0=options.N0, ...
                b=options.b, zDomain=options.zDomain, nModes=options.nModes, ...
                normalization=options.normalization, metadata=options.metadata);
        end
    end

    methods (Static, Access = private)
        function modeNumber = resolveModeNumber(requestedModeNumber, nModes)
            if isempty(requestedModeNumber)
                modeNumber = 1:nModes;
                return;
            end
            modeNumber = reshape(requestedModeNumber,1,[]);
            if length(modeNumber) ~= nModes
                error("IMBasisSet:InvalidModeNumber", ...
                    "modeNumber must have one entry for each retained mode.");
            end
        end

        function normalization = resolveDefaultNormalization(requestedNormalization, evp)
            if ~isempty(requestedNormalization)
                normalization = requestedNormalization;
                return;
            end
            if ~isempty(evp) && isprop(evp, "defaultNormalization") && ~isempty(evp.defaultNormalization)
                normalization = evp.defaultNormalization;
                return;
            end
            normalization = Normalization.unity;
        end
    end
end
