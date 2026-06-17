classdef IMAnalyticalInternalModesBasis
    % Store exact internal-mode functions from an analytical solution family.
    %
    % `IMAnalyticalInternalModesBasis` has the same user-facing `F(z)`,
    % `G(z)`, normalization, and Gram-matrix idiom as numerical internal-mode
    % bases, but it evaluates closed-form functions supplied by an
    % `IMAnalyticalSolution`.
    %
    % ```matlab
    % solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
    % basisSet = solution.internalModes(evp,nModes=4);
    % G = basisSet.G(linspace(-5000,0,128).');
    % ```
    %
    % - Topic: Evaluate analytical modes
    % - Topic: Analyze Gram matrices
    % - Topic: Inspect analytical modes
    % - Topic: Developer topics
    % - Declaration: classdef IMAnalyticalInternalModesBasis

    properties
        % Active normalization rule name.
        %
        % - Topic: Evaluate analytical modes
        normalization = "unity"
    end

    properties (SetAccess = private)
        % Analytical solution family that created this basis.
        %
        % - Topic: Inspect analytical modes
        solution

        % Internal-mode EVP represented by these exact functions.
        %
        % - Topic: Inspect analytical modes
        evp

        % Equivalent depths for retained modes.
        %
        % - Topic: Inspect analytical modes
        h

        % Retained eigenvalues.
        %
        % - Topic: Inspect analytical modes
        eigenvalues

        % Retained mode labels.
        %
        % - Topic: Inspect analytical modes
        modeNumber

        % Physical vertical domain.
        %
        % - Topic: Inspect analytical modes
        zDomain

        % Buoyancy frequency squared function.
        %
        % - Topic: Inspect analytical modes
        N2

        % Additional creation metadata.
        %
        % - Topic: Inspect analytical modes
        metadata

        % Mode-selection diagnostics.
        %
        % Analytical bases do not run a discrete mode-selection pass, so this
        % is usually empty.
        %
        % - Topic: Inspect analytical modes
        modeSelectionDiagnostics
    end

    properties (Access = private)
        rawVariableFunction
        rawUzFunction
        normalizationNameMap
    end

    methods
        function self = IMAnalyticalInternalModesBasis(options)
            % Create an exact internal-mode basis.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: basisSet = IMAnalyticalInternalModesBasis(options)
            % - Parameter options.solution: analytical solution family
            % - Parameter options.evp: internal-mode EVP
            % - Parameter options.h: equivalent depths
            % - Parameter options.modeNumber: retained-mode labels
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.rawVariableFunction: exact `F`/`G` evaluator
            % - Parameter options.rawUzFunction: exact solved-variable derivative evaluator
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: creation metadata
            % - Returns basisSet: exact internal-mode basis
            arguments
                options.solution IMAnalyticalSolution
                options.evp IMInternalModes
                options.h double
                options.modeNumber double
                options.N2 (1,1) function_handle
                options.rawVariableFunction (1,1) function_handle
                options.rawUzFunction (1,1) function_handle
                options.normalization = []
                options.metadata struct = struct()
                options.modeSelectionDiagnostics struct = struct()
            end

            self.solution = options.solution;
            self.evp = options.evp;
            self.h = reshape(options.h,1,[]);
            self.eigenvalues = 1 ./ self.h;
            self.modeNumber = reshape(options.modeNumber,1,[]);
            self.zDomain = options.solution.zDomain;
            self.N2 = options.N2;
            self.metadata = options.metadata;
            self.modeSelectionDiagnostics = options.modeSelectionDiagnostics;
            self.rawVariableFunction = options.rawVariableFunction;
            self.rawUzFunction = options.rawUzFunction;
            self.normalizationNameMap = configureDictionary("string","cell");
            self = self.installNormalizationRules();
            if isempty(options.normalization)
                self.normalization = self.defaultNormalization();
            else
                self.normalization = options.normalization;
            end
        end

        function self = addNormalization(self, name, rule)
            % Add a named normalization rule.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: basisSet = addNormalization(basisSet,name,rule)
            % - Parameter name: normalization rule name
            % - Parameter rule: function handle with signature `scale = rule(basisSet,iMode)`
            % - Returns basisSet: basis set with the rule installed
            arguments
                self IMAnalyticalInternalModesBasis
                name {mustBeTextScalar}
                rule (1,1) function_handle
            end

            name = self.normalizationName(name);
            self.normalizationNameMap{name} = rule;
        end

        function names = normalizationNames(self)
            % Return installed normalization rule names.
            %
            % - Topic: Inspect analytical modes
            % - Declaration: names = normalizationNames(basisSet)
            % - Returns names: installed rule names
            names = sort(string(keys(self.normalizationNameMap)));
        end

        function F = F(self, z, options)
            % Evaluate exact `F` modes.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: F = F(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule
            % - Returns F: exact `F` values
            arguments
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            F = self.rawVariable("F", z) ./ self.normalizationFactors(options.normalization);
        end

        function G = G(self, z, options)
            % Evaluate exact `G` modes.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: G = G(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule
            % - Returns G: exact `G` values
            arguments
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            G = self.rawVariable("G", z) ./ self.normalizationFactors(options.normalization);
        end

        function values = uz(self, z, options)
            % Evaluate exact solved-variable derivatives.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: values = uz(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule
            % - Returns values: exact derivative values
            arguments
                self IMAnalyticalInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            values = self.rawUz(z) ./ self.normalizationFactors(options.normalization);
        end

        function factors = normalizationFactors(self, normalization)
            % Return scale factors for a normalization rule.
            %
            % - Topic: Evaluate analytical modes
            % - Declaration: factors = normalizationFactors(basisSet,normalization)
            % - Parameter normalization: normalization rule name
            % - Returns factors: row vector of scale factors
            arguments
                self IMAnalyticalInternalModesBasis
                normalization = self.normalization
            end

            name = self.normalizationName(normalization);
            if ~isKey(self.normalizationNameMap, name)
                error("IMAnalyticalInternalModesBasis:UnsupportedNormalization", "The analytical basis does not define a ""%s"" normalization.", name);
            end
            rule = self.normalizationNameMap{name};
            factors = zeros(1,length(self.h));
            for iMode = 1:length(self.h)
                factors(iMode) = rule(self, iMode);
            end
            factors = abs(factors);
            factors(factors == 0 | ~isfinite(factors)) = 1;
        end

        function gram = gramMatrix(self, options)
            % Return a Gram matrix for exact `F` or `G` modes.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = gramMatrix(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds
            % - Returns gram: Gram matrix
            arguments
                self IMAnalyticalInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            self.validateZBounds(options.zBounds);
            gram = self.variableGramMatrix(string(options.variable), options.zBounds, true);
        end

        function windowModes = partialWindowModes(self, options)
            % Diagonalize a partial-depth Gram matrix.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: windowModes = partialWindowModes(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds
            % - Returns windowModes: window-mode decomposition
            arguments
                self IMAnalyticalInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            gram = self.gramMatrix(variable=options.variable, zBounds=options.zBounds);
            gram = 0.5*(gram + gram.');
            [R, D] = eig(gram);
            [eigenvalues, sortIndex] = sort(diag(D), "descend");
            windowModes.rotation = R(:,sortIndex);
            windowModes.eigenvalues = eigenvalues(:).';
            windowModes.gramMatrix = gram;
        end

        function spectrum = spectrum(self, coefficients, options)
            % Compute a modal spectrum.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = spectrum(basisSet,coefficients,options)
            % - Parameter coefficients: modal coefficients
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns spectrum: modal spectrum
            arguments
                self IMAnalyticalInternalModesBasis
                coefficients (:,1) double
                options.variable {mustBeTextScalar} = self.evp.formulation
            end

            spectrum = self.crossSpectrum(coefficients, coefficients, variable=options.variable);
        end

        function spectrum = crossSpectrum(self, coefficientsA, coefficientsB, options)
            % Compute a modal cross-spectrum.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB,options)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMAnalyticalInternalModesBasis
                coefficientsA (:,1) double
                coefficientsB (:,1) double
                options.variable {mustBeTextScalar} = self.evp.formulation
            end

            self.validateCoefficientVector(coefficientsA, "coefficientsA");
            self.validateCoefficientVector(coefficientsB, "coefficientsB");
            gram = self.gramMatrix(variable=options.variable);
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
        end
    end

    methods (Hidden)
        function factor = innerProductNormFactor(self, iMode, options)
            % Return the raw inner-product norm factor.
            %
            % - Topic: Developer topics
            % - Declaration: factor = innerProductNormFactor(basisSet,iMode,options)
            % - Parameter iMode: retained mode index
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns factor: normalization factor
            % - Developer: true
            arguments
                self IMAnalyticalInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            gram = self.variableGramMatrix(string(options.variable), self.zDomain, false);
            factor = sqrt(abs(gram(iMode,iMode)));
        end

        function factor = geostrophicNormFactor(self, iMode)
            % Return the geostrophic normalization factor.
            %
            % - Topic: Developer topics
            % - Declaration: factor = geostrophicNormFactor(basisSet,iMode)
            % - Returns factor: normalization factor
            % - Developer: true
            if self.modeNumber(iMode) == 0
                gram = self.variableGramMatrix("F", self.zDomain, false);
                factor = sqrt(abs(gram(iMode,iMode)))/sqrt(diff(self.zDomain));
                return;
            end
            factor = self.innerProductNormFactor(iMode, variable="G");
        end

        function factor = maxAbsFactor(self, iMode, options)
            % Return the maximum amplitude of `F` or `G`.
            %
            % - Topic: Developer topics
            % - Declaration: factor = maxAbsFactor(basisSet,iMode,options)
            % - Parameter iMode: retained mode index
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns factor: maximum amplitude
            % - Developer: true
            arguments
                self IMAnalyticalInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            values = self.rawVariable(string(options.variable), self.integrationGrid(self.zDomain));
            factor = max(abs(values(:,iMode)));
        end

        function factor = surfacePressureNormFactor(self, iMode)
            % Return the raw surface `F` value.
            %
            % - Topic: Developer topics
            % - Declaration: factor = surfacePressureNormFactor(basisSet,iMode)
            % - Returns factor: normalization factor
            % - Developer: true
            values = self.rawVariable("F", self.zDomain(2));
            factor = values(1,iMode);
            if ~isfinite(factor) || abs(factor) <= 1e-12
                factor = 1;
            end
        end
    end

    methods (Access = private)
        function self = installNormalizationRules(self)
            self = self.addNormalization("unity", @(basisSet,iMode) basisSet.innerProductNormFactor(iMode));
            self = self.addNormalization("uMax", @(basisSet,iMode) basisSet.maxAbsFactor(iMode, variable="F"));
            self = self.addNormalization("wMax", @(basisSet,iMode) basisSet.maxAbsFactor(iMode, variable="G"));
            self = self.addNormalization("surfacePressure", @(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode));
            if self.evp.modeFamily == "geostrophic"
                self = self.addNormalization("geostrophic", @(basisSet,iMode) basisSet.geostrophicNormFactor(iMode));
            end
            if string(self.evp.name) == "waveModesAtWavenumber"
                self = self.addNormalization("kConstant", @(basisSet,iMode) basisSet.innerProductNormFactor(iMode, variable="G"));
            end
        end

        function name = defaultNormalization(self)
            if self.evp.modeFamily == "geostrophic"
                name = "geostrophic";
            elseif string(self.evp.name) == "waveModesAtWavenumber"
                name = "kConstant";
            else
                name = "unity";
            end
        end

        function values = rawVariable(self, variable, z)
            values = self.rawVariableFunction(self, string(variable), z(:));
        end

        function values = rawUz(self, z)
            values = self.rawUzFunction(self, z(:));
        end

        function gram = variableGramMatrix(self, variable, zBounds, useNormalized)
            spec = self.evp.innerProduct(variable);
            IMAnalyticalInternalModesBasis.assertInnerProductAvailable(spec);
            z = self.integrationGrid(zBounds);
            if useNormalized
                values = self.rawVariable(variable, z) ./ self.normalizationFactors(self.normalization);
            else
                values = self.rawVariable(variable, z);
            end
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, self.context());
            if isscalar(weight)
                weight = weight*ones(size(z));
            end
            gram = zeros(size(values,2), size(values,2));
            for iMode = 1:size(values,2)
                for jMode = iMode:size(values,2)
                    value = trapz(z, weight(:).*values(:,iMode).*values(:,jMode));
                    if variable == self.evp.formulation
                        endpointWeights = [spec.surfaceWeights; spec.bottomWeights];
                        for iWeight = 1:numel(endpointWeights)
                            value = value + self.endpointWeightContribution(endpointWeights(iWeight), iMode, jMode, useNormalized, zBounds);
                        end
                    end
                    for iTerm = 1:numel(spec.endpointInnerProductTerms)
                        value = value + self.endpointTermContribution(spec.endpointInnerProductTerms(iTerm), iMode, jMode, useNormalized, zBounds);
                    end
                    gram(iMode,jMode) = value;
                    gram(jMode,iMode) = value;
                end
            end
        end

        function value = endpointWeightContribution(self, weight, iMode, jMode, useNormalized, zBounds)
            zEndpoint = self.endpointZ(weight.location);
            if ~self.boundsIncludeEndpoint(zBounds, zEndpoint)
                value = 0;
                return;
            end
            if useNormalized
                uValues = self.rawVariable(self.evp.formulation, zEndpoint) ./ self.normalizationFactors(self.normalization);
                uzValues = self.rawUz(zEndpoint) ./ self.normalizationFactors(self.normalization);
            else
                uValues = self.rawVariable(self.evp.formulation, zEndpoint);
                uzValues = self.rawUz(zEndpoint);
            end
            pValue = IMEigenvalueProblem.evaluateCoefficient(self.evp.p, zEndpoint, self.context());
            values = weight.c*uValues - weight.d*pValue*uzValues;
            value = weight.coefficient*values(iMode)*values(jMode);
        end

        function value = endpointTermContribution(self, term, iMode, jMode, useNormalized, zBounds)
            zEndpoint = self.endpointZ(term.location);
            if ~self.boundsIncludeEndpoint(zBounds, zEndpoint)
                value = 0;
                return;
            end
            if useNormalized
                values = self.rawVariable(term.variable, zEndpoint) ./ self.normalizationFactors(self.normalization);
            else
                values = self.rawVariable(term.variable, zEndpoint);
            end
            value = term.coefficient*values(iMode)*values(jMode);
        end

        function z = integrationGrid(~, zBounds)
            z = linspace(min(zBounds), max(zBounds), 1024).';
        end

        function context = context(self)
            context.zDomain = self.zDomain;
            context.coordinateKind = "analyticalSolution";
            context.N2 = @(z) self.N2(z);
            context.f0 = self.evp.f0;
            context.g = self.evp.g;
            context.formulation = self.evp.formulation;
            parameterFields = fieldnames(self.evp.parameters);
            for iField = 1:numel(parameterFields)
                fieldName = parameterFields{iField};
                context.(fieldName) = self.evp.parameters.(fieldName);
            end
        end

        function validateZBounds(~, zBounds)
            if zBounds(1) >= zBounds(2)
                error("IMAnalyticalInternalModesBasis:InvalidInterval", "zBounds must be increasing.");
            end
        end

        function validateCoefficientVector(self, coefficients, argumentName)
            if any(~isfinite(coefficients(:)))
                error("IMAnalyticalInternalModesBasis:InvalidCoefficients", "%s must contain finite values.", char(argumentName));
            end
            if length(coefficients) ~= length(self.h)
                error("IMAnalyticalInternalModesBasis:InvalidCoefficientCount", "%s must contain one value for each retained mode (%d).", char(argumentName), length(self.h));
            end
        end

        function zEndpoint = endpointZ(self, location)
            switch string(location)
                case "surface"
                    zEndpoint = self.zDomain(2);
                case "bottom"
                    zEndpoint = self.zDomain(1);
                otherwise
                    error("IMAnalyticalInternalModesBasis:InvalidBoundaryLocation", "Boundary location must be ""surface"" or ""bottom"".");
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
    end

    methods (Static, Access = private)
        function assertInnerProductAvailable(spec)
            if isfield(spec, "hasInnerProduct") && spec.hasInnerProduct
                return;
            end
            error("IMAnalyticalInternalModesBasis:UnavailableInnerProduct", "The %s inner product is unavailable for this analytical basis. %s", string(spec.variable), string(spec.reason));
        end
    end
end
