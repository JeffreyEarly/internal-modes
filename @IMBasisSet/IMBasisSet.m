classdef IMBasisSet
    % Store solved scalar canonical EVP modes.
    %
    % `IMBasisSet` stores the scalar modes selected from a canonical EVP and
    % evaluates the solved variable `u` and derivative `uz`. Modal
    % normalization is applied lazily when values, Gram matrices, or spectra
    % are requested. For each retained mode,
    % $$u_j^{\mathrm{out}}(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
    % where the scale factor $$s_j$$ is supplied by the active
    % normalization rule.
    %
    % ```matlab
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet.normalization = "unity";
    % u = basisSet.u(z);
    % factors = basisSet.normalizationFactors("unity");
    % ```
    %
    % - Topic: Create basis sets
    % - Topic: Evaluate basis sets
    % - Topic: Analyze Gram matrices
    % - Topic: Inspect basis sets
    % - Topic: Developer topics
    % - Declaration: classdef IMBasisSet

    properties
        % Active normalization rule name.
        %
        % This string selects a field in `basisSet.evp.normalizationRules`;
        % create custom rules on the EVP, not on the basis set.
        % The selected rule returns the scale factor $$s_j$$ used by `u`,
        % `uz`, and Gram-matrix methods. Passing `normalization=...` to an
        % evaluation method overrides this property for that call.
        %
        % - Topic: Evaluate basis sets
        normalization = "unity"
    end

    properties (SetAccess = protected)
        % Solver that created the native modes.
        %
        % This is the discretization object used to interpolate native
        % coefficient columns. It is mainly useful for diagnostics and
        % developer workflows.
        %
        % - Topic: Developer topics
        % - Developer: true
        solver

        % EVP descriptor that was solved.
        %
        % - Topic: Inspect basis sets
        evp

        % Native mode columns.
        %
        % These are the columns returned by the numerical solver before
        % interpolation onto physical coordinates and before modal
        % normalization. Most workflows should call `u(z)` or `uz(z)`
        % instead.
        %
        % - Topic: Developer topics
        % - Developer: true
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

        % Mode-selection diagnostics.
        %
        % This is the diagnostics struct returned by
        % `evp.modeSelectionDiagnostics(...)` when the solver selected and
        % labeled retained modes. Numerical solves use it to record
        % negative-mode bounds and zero-mode status; analytical basis sets
        % may leave it as an empty struct.
        %
        % - Topic: Inspect basis sets
        modeSelectionDiagnostics

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
            % - Parameter options.modeSelectionDiagnostics: mode-selection diagnostics
            % - Parameter options.normalization: active normalization rule name
            % - Parameter options.metadata: additional metadata
            % - Parameter options.zDomain: physical vertical domain
            % - Returns basisSet: solved scalar basis set
            arguments
                options.solver = []
                options.evp = []
                options.nativeModes double = zeros(0,0)
                options.eigenvalues double = zeros(1,0)
                options.modeNumber double = []
                options.modeSelectionDiagnostics struct = struct()
                options.normalization = []
                options.metadata struct = struct()
                options.zDomain (1,2) double = [NaN NaN]
            end

            self.solver = options.solver;
            self.evp = options.evp;
            self.nativeModes = options.nativeModes;
            self.eigenvalues = reshape(options.eigenvalues,1,[]);
            nModes = max(size(options.nativeModes,2), length(self.eigenvalues));
            self.modeNumber = IMBasisSet.resolveModeNumber(options.modeNumber, nModes);
            self.modeSelectionDiagnostics = options.modeSelectionDiagnostics;
            self.normalization = IMBasisSet.resolveDefaultNormalization(options.normalization, options.evp);
            self.metadata = options.metadata;
            self.zDomain = options.zDomain;

            if any(isnan(self.zDomain)) && ~isempty(self.evp)
                self.zDomain = self.evp.zDomain;
            elseif ~isempty(self.solver)
                if any(isnan(self.zDomain))
                    self.zDomain = self.solver.zDomain;
                end
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
                z (:,1) double {mustBeReal, mustBeFinite}
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
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            values = self.rawUz(z) ./ self.normalizationFactors(options.normalization);
        end

        function factors = normalizationFactors(self, normalization)
            % Return scale factors for a normalization rule.
            %
            % For a requested rule name this returns the row vector
            % $$s_j$$ used by evaluation methods:
            % $$u_j^{\mathrm{out}}=u_j^{\mathrm{raw}}/s_j.$$
            % The default scalar `"unity"` rule uses
            % $$s_j=\sqrt{|\langle u_j,u_j\rangle|}$$ with the EVP inner
            % product. Custom rules are created on
            % `basisSet.evp.normalizationRules`, the default rule is named
            % by `basisSet.evp.defaultNormalization`, and the active rule
            % is selected by `basisSet.normalization`.
            %
            % ```matlab
            % factors = basisSet.normalizationFactors("unity");
            % uUnity = basisSet.u(z,normalization="unity");
            % ```
            %
            % - Topic: Evaluate basis sets
            % - Declaration: factors = normalizationFactors(basisSet,normalization)
            % - Parameter normalization: normalization convention
            % - Returns factors: row vector of positive scale factors
            arguments
                self IMBasisSet
                normalization = self.normalization
            end

            name = char(self.normalizationName(normalization));
            if isempty(self.evp) || ~isfield(self.evp.normalizationRules, name)
                error("IMBasisSet:UnsupportedNormalization", ...
                    "The EVP does not define a ""%s"" normalization.", name);
            end
            normalizeMode = self.evp.normalizationRules.(name);
            nModes = self.retainedModeCount();
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
            % - Topic: Developer topics
            % - Declaration: modes = normalizedNativeModes(basisSet,normalization)
            % - Parameter normalization: normalization convention
            % - Returns modes: scaled native mode columns
            % - Developer: true
            arguments
                self IMBasisSet
                normalization = self.normalization
            end

            factors = self.normalizationFactors(normalization);
            modes = self.nativeModes ./ factors;
        end

        function gram = gramMatrix(self)
            % Return the full-domain scalar Gram matrix.
            %
            % For normalized scalar modes, entries are
            % $$M_{ij}=\int r u_i u_j\,dz+
            % \sum_\ell \gamma_\ell L_\ell[u_i]L_\ell[u_j].$$
            % The endpoint terms are the active metric weights from
            % `evp.innerProduct()`.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = gramMatrix(basisSet)
            % - Returns gram: scalar Gram matrix
            arguments
                self IMBasisSet
            end

            gram = self.partialGramMatrix(self.zDomain(1), self.zDomain(2));
        end

        function gram = partialGramMatrix(self, zMin, zMax)
            % Return a partial-domain scalar Gram matrix.
            %
            % Interior integrals are restricted to `[zMin,zMax]`. Endpoint
            % metric terms are included only when the requested interval
            % contains the corresponding physical endpoint.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = partialGramMatrix(basisSet,zMin,zMax)
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns gram: scalar Gram matrix
            arguments
                self IMBasisSet
                zMin (1,1) double {mustBeReal, mustBeFinite}
                zMax (1,1) double {mustBeReal, mustBeFinite}
            end

            self.validateZBounds(zMin, zMax);
            gram = self.scalarGramMatrix([zMin zMax], true);
        end

        function windowModes = partialWindowModes(self, zMin, zMax)
            % Diagonalize a partial scalar Gram matrix.
            %
            % This computes the eigendecomposition of the symmetric
            % partial-domain Gram matrix and sorts window-mode eigenvalues
            % from largest to smallest.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: windowModes = partialWindowModes(basisSet,zMin,zMax)
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns windowModes: window-mode decomposition
            arguments
                self IMBasisSet
                zMin (1,1) double {mustBeReal, mustBeFinite}
                zMax (1,1) double {mustBeReal, mustBeFinite}
            end

            gram = self.partialGramMatrix(zMin, zMax);
            gram = 0.5*(gram + gram.');
            [R, D] = eig(gram);
            [eigenvalues, sortIndex] = sort(diag(D), "descend");
            windowModes.rotation = R(:,sortIndex);
            windowModes.eigenvalues = eigenvalues(:).';
            windowModes.gramMatrix = gram;
        end

        function spectrum = spectrum(self, coefficients)
            % Compute a scalar modal spectrum.
            %
            % For modal coefficients $$c_j$$ this returns
            % $$S_j=M_{jj}|c_j|^2,$$ where $$M$$ is the full-domain scalar
            % Gram matrix.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = spectrum(basisSet,coefficients)
            % - Parameter coefficients: modal coefficients
            % - Returns spectrum: modal spectrum
            arguments
                self IMBasisSet
                coefficients (:,1) double
            end

            self.validateCoefficientVector(coefficients, "coefficients");
            spectrum = self.crossSpectrum(coefficients, coefficients);
        end

        function spectrum = crossSpectrum(self, coefficientsA, coefficientsB)
            % Compute a scalar modal cross-spectrum.
            %
            % For modal coefficient vectors $$a_j$$ and $$b_j$$ this
            % returns $$S_j=M_{jj}\operatorname{Re}(a_j b_j^*)$$.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMBasisSet
                coefficientsA (:,1) double
                coefficientsB (:,1) double
            end

            self.validateCoefficientVector(coefficientsA, "coefficientsA");
            self.validateCoefficientVector(coefficientsB, "coefficientsB");
            gram = self.gramMatrix();
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
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

        function factor = innerProductNormFactor(self, iMode)
            % Return the scalar inner-product norm factor.
            %
            % This is the raw factor
            % $$s_j=\sqrt{|\langle u_j,u_j\rangle|}$$ computed before any
            % modal normalization has been applied.
            %
            % - Topic: Developer topics
            % - Declaration: factor = innerProductNormFactor(basisSet,iMode)
            % - Developer: true
            arguments
                self IMBasisSet
                iMode (1,1) double {mustBeInteger, mustBePositive}
            end

            gram = self.scalarGramMatrix(self.zDomain, false);
            factor = sqrt(abs(gram(iMode,iMode)));
        end

        function factor = maxAbsFactor(self, iMode)
            % Return the maximum scalar amplitude.
            %
            % This is $$s_j=\max_z |u_j^{\mathrm{raw}}(z)|$$ on the
            % basis-set integration grid.
            %
            % - Topic: Developer topics
            % - Declaration: factor = maxAbsFactor(basisSet,iMode)
            % - Developer: true
            arguments
                self IMBasisSet
                iMode (1,1) double {mustBeInteger, mustBePositive}
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
                nGrid = max(256, 4*max(1,self.retainedModeCount()));
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
            spec = self.evp.innerProduct();
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

        function nModes = retainedModeCount(self)
            nModes = max([length(self.eigenvalues), size(self.nativeModes,2), length(self.modeNumber)]);
        end

        function validateCoefficientVector(self, coefficients, argumentName)
            if any(~isfinite(coefficients(:)))
                error("IMBasisSet:InvalidCoefficients", ...
                    "%s must contain finite values.", char(argumentName));
            end
            nModes = self.retainedModeCount();
            if length(coefficients) ~= nModes
                error("IMBasisSet:InvalidCoefficientCount", ...
                    "%s must contain one value for each retained mode (%d).", char(argumentName), nModes);
            end
        end

        function validateZBounds(~, zMin, zMax)
            if zMin >= zMax
                error("IMBasisSet:InvalidInterval", ...
                    "zMin must be less than zMax.");
            end
        end

        function isIncluded = boundsIncludeEndpoint(~, zBounds, endpoint)
            tolerance = 100*eps(max(1,max(abs([zBounds(:); endpoint]))));
            isIncluded = abs(min(zBounds) - endpoint) <= tolerance || abs(max(zBounds) - endpoint) <= tolerance;
        end

        function name = normalizationName(~, normalization)
            parts = split(string(normalization), ".");
            name = parts(end);
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
            normalization = "unity";
        end
    end
end
