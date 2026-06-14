classdef IMInternalModesBasis < IMBasisSet
    % Store solved internal-mode basis functions.
    %
    % `IMInternalModesBasis` evaluates the physical `F` and `G` variables
    % from the solved canonical scalar mode. If the EVP solves `G`, then
    % $$F_j=h_j G'_j.$$ If the EVP solves `F`, then
    % $$G_j=-gN^{-2}F'_j.$$
    % Normalization is shared across both variables: if a rule gives scale
    % $$s_j$$, then both diagnostic variables for mode $$j$$ are divided by
    % that same factor.
    %
    % ```matlab
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet.normalization = Normalization.geostrophic;
    % F = basisSet.F(z);
    % G = basisSet.G(z);
    % ```
    %
    % - Topic: Evaluate internal-mode bases
    % - Topic: Analyze Gram matrices
    % - Topic: Inspect internal-mode bases
    % - Topic: Developer topics
    % - Declaration: classdef IMInternalModesBasis < IMBasisSet

    properties (SetAccess = protected)
        % Equivalent depths for the retained internal modes.
        %
        % For numerical internal-mode solves, these are computed from the
        % parent EVP as
        % $$h_j=\texttt{evp.hFromEigenvalue}(\lambda_j),$$
        % where $$\lambda_j$$ is `eigenvalues(j)`. Analytical basis classes
        % may pass exact equivalent depths directly.
        %
        % - Topic: Inspect internal-mode bases
        h
    end

    methods
        function self = IMInternalModesBasis(options)
            % Create an internal-mode basis set.
            %
            % - Topic: Evaluate internal-mode bases
            % - Declaration: basisSet = IMInternalModesBasis(options)
            % - Parameter options.solver: solver reference
            % - Parameter options.evp: internal-mode EVP descriptor
            % - Parameter options.nativeModes: native mode columns
            % - Parameter options.eigenvalues: retained eigenvalues
            % - Parameter options.h: equivalent depths
            % - Parameter options.modeNumber: physical mode numbers
            % - Parameter options.index: selection diagnostics
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.N2Function: buoyancy frequency squared function
            % - Returns basisSet: internal-mode basis set
            arguments
                options.solver = []
                options.evp IMInternalModes
                options.nativeModes double = zeros(0,0)
                options.eigenvalues double = zeros(1,0)
                options.h double = zeros(1,0)
                options.modeNumber double = []
                options.index struct = struct()
                options.normalization = []
                options.metadata struct = struct()
                options.zDomain (1,2) double = [NaN NaN]
                options.N2Function = []
            end

            self@IMBasisSet(solver=options.solver, evp=options.evp, ...
                nativeModes=options.nativeModes, eigenvalues=options.eigenvalues, ...
                modeNumber=options.modeNumber, index=options.index, normalization=options.normalization, ...
                metadata=options.metadata, zDomain=options.zDomain, N2Function=options.N2Function);
            self.h = IMInternalModesBasis.resolveEquivalentDepths(options.h, self.eigenvalues, options.evp);
            self.validateEquivalentDepthCount();
        end

        function G = G(self, z, options)
            % Evaluate `G` modes.
            %
            % If the EVP formulation is `G`, this evaluates the solved
            % canonical variable. If the formulation is `F`, `G` is recovered
            % from $$G=-gN^{-2}F_z$$.
            %
            % - Topic: Evaluate internal-mode bases
            % - Declaration: G = G(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization to apply
            % - Returns G: evaluated `G` modes
            arguments
                self IMInternalModesBasis
                z (:,1) double
                options.normalization = self.normalization
            end

            G = self.rawVariable("G", z) ./ self.normalizationFactors(options.normalization);
        end

        function F = F(self, z, options)
            % Evaluate `F` modes.
            %
            % If the EVP formulation is `F`, this evaluates the solved
            % canonical variable. If the formulation is `G`, `F` is recovered
            % from $$F=hG_z$$.
            %
            % - Topic: Evaluate internal-mode bases
            % - Declaration: F = F(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization to apply
            % - Returns F: evaluated `F` modes
            arguments
                self IMInternalModesBasis
                z (:,1) double
                options.normalization = self.normalization
            end

            F = self.rawVariable("F", z) ./ self.normalizationFactors(options.normalization);
        end

        function values = evaluate(self, variable, z, options)
            % Evaluate `F` or `G`.
            %
            % This is the variable-dispatch form of `F(z)` and `G(z)`.
            %
            % - Topic: Evaluate internal-mode bases
            % - Declaration: values = evaluate(basisSet,variable,z,options)
            % - Parameter variable: `"F"` or `"G"`
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization to apply
            % - Returns values: evaluated modes
            arguments
                self IMInternalModesBasis
                variable {mustBeTextScalar}
                z (:,1) double
                options.normalization = self.normalization
            end

            switch IMInternalModesBasis.validateVariable(variable)
                case "G"
                    values = self.G(z, normalization=options.normalization);
                case "F"
                    values = self.F(z, normalization=options.normalization);
            end
        end

        function values = evaluateAll(self, z)
            % Evaluate both `F` and `G`.
            %
            % - Topic: Evaluate internal-mode bases
            % - Declaration: values = evaluateAll(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: structure with fields `F` and `G`
            values.F = self.F(z);
            values.G = self.G(z);
        end

        function gram = gramMatrix(self, variable)
            % Return the full-depth Gram matrix for `F` or `G`.
            %
            % The matrix uses `evp.innerProduct(variable)`. For `G`, the
            % interior weight is $$N^2/g$$; for `F`, it is one. Endpoint
            % terms are included when the requested variable is the solved
            % canonical formulation and the endpoint condition supplies
            % metric weights.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = gramMatrix(basisSet,variable)
            % - Parameter variable: `"F"` or `"G"`
            % - Returns gram: Gram matrix
            arguments
                self IMInternalModesBasis
                variable {mustBeTextScalar} = self.evp.formulation
            end

            gram = self.partialGramMatrix(variable, self.zDomain(1), self.zDomain(2));
        end

        function gram = partialGramMatrix(self, varargin)
            % Return a partial-depth Gram matrix for `F` or `G`.
            %
            % Interior integrals are restricted to `[zMin,zMax]`; endpoint
            % metric terms are included only when the interval contains the
            % corresponding endpoint. If `variable` is omitted, the solved
            % formulation is used.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = partialGramMatrix(basisSet,variable,zMin,zMax)
            % - Parameter variable: optional variable name, `"F"` or `"G"`
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns gram: Gram matrix
            [variable, zMin, zMax] = self.parseVariableBounds(varargin{:});
            gram = self.variableGramMatrix(variable, [zMin zMax], true);
        end

        function windowModes = partialWindowModes(self, varargin)
            % Diagonalize a partial-depth Gram matrix for `F` or `G`.
            %
            % If `variable` is omitted, the solved formulation is used.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: windowModes = partialWindowModes(basisSet,variable,zMin,zMax)
            % - Parameter variable: optional variable name, `"F"` or `"G"`
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns windowModes: window-mode decomposition
            [variable, zMin, zMax] = self.parseVariableBounds(varargin{:});
            gram = self.variableGramMatrix(variable, [zMin zMax], true);
            gram = 0.5*(gram + gram.');
            [R, D] = eig(gram);
            [eigenvalues, sortIndex] = sort(diag(D), "descend");
            windowModes.rotation = R(:,sortIndex);
            windowModes.eigenvalues = eigenvalues(:).';
            windowModes.gramMatrix = gram;
        end

        function spectrum = spectrum(self, coefficients, options)
            % Compute an internal-mode modal spectrum.
            %
            % If `options.variable` is omitted, the solved formulation is
            % used.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = spectrum(basisSet,coefficients,options)
            % - Parameter coefficients: modal coefficients
            % - Parameter options.variable: optional variable name, `"F"` or `"G"`
            % - Returns spectrum: modal spectrum
            arguments
                self IMInternalModesBasis
                coefficients (:,1) double
                options.variable {mustBeTextScalar} = self.evp.formulation
            end

            spectrum = self.crossSpectrum(coefficients, coefficients, variable=options.variable);
        end

        function spectrum = crossSpectrum(self, coefficientsA, coefficientsB, options)
            % Compute an internal-mode modal cross-spectrum.
            %
            % If `options.variable` is omitted, the solved formulation is
            % used.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB,options)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Parameter options.variable: optional variable name, `"F"` or `"G"`
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMInternalModesBasis
                coefficientsA (:,1) double
                coefficientsB (:,1) double
                options.variable {mustBeTextScalar} = self.evp.formulation
            end

            gram = self.gramMatrix(options.variable);
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
        end
    end

    methods (Hidden)
        function self = orientModeSigns(self)
            % Orient modes so the surface `F` value is positive when possible.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = orientModeSigns(basisSet)
            % - Developer: true
            if isempty(self.nativeModes)
                return;
            end
            zSurface = self.zDomain(2);
            zGrid = self.innerProductGrid(self.zDomain);
            FSurface = self.rawVariableForSign("F", zSurface);
            FGrid = self.rawVariableForSign("F", zGrid);
            GSurface = self.rawVariableForSign("G", zSurface);
            GGrid = self.rawVariableForSign("G", zGrid);
            signs = ones(1,size(self.nativeModes,2));
            for iMode = 1:size(self.nativeModes,2)
                reference = self.referenceValue(FSurface, FGrid, iMode);
                if ~isfinite(reference) || reference == 0
                    reference = self.referenceValue(GSurface, GGrid, iMode);
                end
                if isfinite(reference) && reference < 0
                    signs(iMode) = -1;
                end
            end
            self.nativeModes = self.nativeModes .* signs;
        end

        function factor = innerProductNormFactor(self, variable, iMode)
            % Return the `F` or `G` inner-product norm factor.
            %
            % This is the raw factor
            % $$s_j=\sqrt{|\langle V_j,V_j\rangle|}$$ for `variable` equal
            % to `F` or `G`.
            %
            % - Topic: Developer topics
            % - Declaration: factor = innerProductNormFactor(basisSet,variable,iMode)
            % - Developer: true
            variable = IMInternalModesBasis.validateVariable(variable);
            gram = self.variableGramMatrix(variable, self.zDomain, false);
            factor = sqrt(abs(gram(iMode,iMode)));
        end

        function factor = geostrophicNormFactor(self, iMode)
            % Return the hydrostatic geostrophic normalization factor.
            %
            % Baroclinic modes use the `G` inner-product norm. A barotropic
            % zero mode uses the `F` norm divided by
            % $$\sqrt{z_\mathrm{surface}-z_\mathrm{bottom}}$$.
            %
            % - Topic: Developer topics
            % - Declaration: factor = geostrophicNormFactor(basisSet,iMode)
            % - Developer: true
            if self.modeNumber(iMode) == 0
                gram = self.variableGramMatrix("F", self.zDomain, false);
                factor = sqrt(abs(gram(iMode,iMode)))/sqrt(diff(self.zDomain));
                return;
            end
            factor = self.innerProductNormFactor("G", iMode);
        end

        function factor = maxAbsFactor(self, variable, iMode)
            % Return the maximum amplitude of `F` or `G`.
            %
            % This is $$s_j=\max_z |V_j^{\mathrm{raw}}(z)|$$ for the
            % requested variable on the basis-set integration grid.
            %
            % - Topic: Developer topics
            % - Declaration: factor = maxAbsFactor(basisSet,variable,iMode)
            % - Developer: true
            variable = IMInternalModesBasis.validateVariable(variable);
            z = self.integrationGrid(self.zDomain);
            values = self.rawVariable(variable, z);
            factor = max(abs(values(:,iMode)));
        end

        function factor = surfacePressureNormFactor(self, iMode)
            % Return the raw surface `F` value.
            %
            % This scale gives unit surface `F` value when the raw surface
            % value is finite and nonzero.
            %
            % - Topic: Developer topics
            % - Declaration: factor = surfacePressureNormFactor(basisSet,iMode)
            % - Developer: true
            values = self.rawVariable("F", self.zDomain(2));
            factor = values(1,iMode);
            if ~isfinite(factor) || abs(factor) <= 1e-12
                factor = 1;
            end
        end
    end

    methods (Access = protected)
        function values = rawU(self, z)
            if ~isempty(self.solver)
                values = self.solver.evaluateNativeModes(self.nativeModes, z);
                return;
            end
            values = self.rawVariable(self.evp.formulation, z);
        end

        function values = rawUz(self, z)
            if ~isempty(self.solver)
                values = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
                return;
            end
            z = z(:);
            if length(z) == 1
                scale = max(1,abs(z));
                dz = sqrt(eps)*scale;
                values = (self.rawU(z + dz) - self.rawU(z - dz))/(2*dz);
            else
                values = gradient(self.rawU(z), z);
            end
        end

        function values = rawVariable(self, variable, z)
            variable = IMInternalModesBasis.validateVariable(variable);
            if variable == self.evp.formulation
                if isempty(self.solver)
                    self.unsupported("evaluate analytical " + variable + " modes without a subclass implementation");
                end
                values = self.solver.evaluateNativeModes(self.nativeModes, z);
                return;
            end

            switch self.evp.formulation
                case "G"
                    if variable ~= "F"
                        self.unsupported("evaluate " + variable);
                    end
                    values = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1) .* self.h;
                case "F"
                    if variable ~= "G"
                        self.unsupported("evaluate " + variable);
                    end
                    dFdz = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
                    values = -(self.evp.g./self.N2(z(:))).*dFdz;
                otherwise
                    self.unsupported("evaluate " + variable);
            end
        end

        function gram = variableGramMatrix(self, variable, zBounds, useNormalized)
            z = self.innerProductGrid(zBounds);
            if useNormalized
                values = self.evaluate(variable, z);
            else
                values = self.rawVariable(variable, z);
            end
            spec = self.evp.innerProduct(variable);
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, self.context());
            if isscalar(weight)
                weight = weight*ones(size(z));
            end
            gram = zeros(size(values,2), size(values,2));
            for iMode = 1:size(values,2)
                for jMode = iMode:size(values,2)
                    integrand = weight(:).*values(:,iMode).*values(:,jMode);
                    value = self.integrateInnerProduct(z, integrand, zBounds);
                    if variable == self.evp.formulation
                        value = value + self.endpointMetricValue([spec.surfaceWeights; spec.bottomWeights], ...
                            iMode, jMode, useNormalized, zBounds);
                    end
                    gram(iMode,jMode) = value;
                    gram(jMode,iMode) = value;
                end
            end
        end
    end

    methods (Access = private)
        function [variable, zMin, zMax] = parseVariableBounds(self, varargin)
            switch numel(varargin)
                case 2
                    variable = self.evp.formulation;
                    zMin = varargin{1};
                    zMax = varargin{2};
                case 3
                    variable = varargin{1};
                    zMin = varargin{2};
                    zMax = varargin{3};
                otherwise
                    error("IMInternalModesBasis:InvalidInput", ...
                        "Use partial Gram-matrix methods as method(zMin,zMax) or method(variable,zMin,zMax).");
            end
            variable = IMInternalModesBasis.validateVariable(variable);
            if ~isnumeric(zMin) || ~isscalar(zMin) || ~isfinite(zMin) || ~isnumeric(zMax) || ~isscalar(zMax) || ~isfinite(zMax)
                error("IMInternalModesBasis:InvalidBounds", ...
                    "zMin and zMax must be finite numeric scalars.");
            end
            zMin = double(zMin);
            zMax = double(zMax);
        end

        function values = rawVariableForSign(self, variable, z)
            values = [];
            try
                values = self.rawVariable(variable, z);
            catch exception
                if string(exception.identifier) ~= "IMBasisSet:UnsupportedOperation"
                    rethrow(exception)
                end
            end
        end

        function value = referenceValue(~, surfaceValues, gridValues, iMode)
            value = NaN;
            tolerance = 0;
            if ~isempty(surfaceValues)
                value = surfaceValues(1,iMode);
                tolerance = 1e-10*max(1,abs(value));
            end
            if (isempty(surfaceValues) || abs(value) <= tolerance) && ~isempty(gridValues)
                [scale, index] = max(abs(gridValues(:,iMode)));
                value = gridValues(index,iMode);
                tolerance = 1e-10*max(1,scale);
            end
            if abs(value) <= tolerance
                value = 0;
            end
        end
    end

    methods (Static, Access = private)
        function h = resolveEquivalentDepths(requestedH, eigenvalues, evp)
            h = reshape(requestedH,1,[]);
            if isempty(h) && ~isempty(eigenvalues)
                h = evp.hFromEigenvalue(reshape(eigenvalues,1,[]));
                h = reshape(h,1,[]);
            end
        end

        function variable = validateVariable(variable)
            variable = string(variable);
            if variable ~= "F" && variable ~= "G"
                error("IMInternalModesBasis:InvalidVariable", ...
                    "variable must be ""F"" or ""G"".");
            end
        end
    end

    methods (Access = private)
        function validateEquivalentDepthCount(self)
            if isempty(self.h)
                return;
            end
            nModes = max([size(self.nativeModes,2), length(self.eigenvalues), length(self.modeNumber)]);
            if length(self.h) ~= nModes
                error("IMInternalModesBasis:InvalidEquivalentDepthCount", ...
                    "h must contain one equivalent depth for each retained mode.");
            end
        end
    end
end
