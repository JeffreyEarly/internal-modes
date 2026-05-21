classdef IMBasisSet
    % Store solved EVP modes with basis-owned normalization.
    %
    % `IMBasisSet` stores the native modes returned by a solver,
    % or exact analytical solution metadata supplied by a subclass.
    % Normalization is applied lazily when modes, spectra, or Gram matrices
    % are requested.
    %
    % ```matlab
    % basisSet.normalization = Normalization.geostrophic;
    % G = basisSet.G(linspace(-1000,0,64).');
    % ```
    %
    % - Topic: Create basis sets
    % - Topic: Evaluate basis sets
    % - Topic: Inspect basis sets
    % - Topic: Analyze Gram matrices
    % - Topic: Developer topics
    % - Declaration: classdef IMBasisSet

    properties
        % Active modal normalization.
        %
        % Changing this property rescales evaluated modes without solving a
        % new EVP. When omitted at construction, this is initialized from
        % the EVP's declared default normalization.
        %
        % - Topic: Evaluate basis sets
        normalization = Normalization.kConstant
    end

    properties (SetAccess = private)
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
        % - Topic: Inspect basis sets
        eigenvalues

        % Equivalent depths.
        %
        % - Topic: Inspect basis sets
        h

        % Physical mode numbers.
        %
        % `-1` identifies a surface boundary mode, `-2` identifies a bottom
        % boundary mode, `0` identifies a true null mode with
        % $$F_0(z)=1$$ and $$G_0(z)=0$$, and positive labels identify
        % interior baroclinic modes.
        %
        % - Topic: Inspect basis sets
        modeNumber

        % Observed and expected eigenvalue index counts.
        %
        % - Topic: Inspect basis sets
        index

        % Additional metadata.
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
            % Create a solved basis set.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet(options)
            % - Parameter options.solver: solver reference
            % - Parameter options.evp: EVP descriptor
            % - Parameter options.nativeModes: native mode columns
            % - Parameter options.eigenvalues: retained eigenvalues
            % - Parameter options.h: equivalent depths
            % - Parameter options.modeNumber: physical mode numbers
            % - Parameter options.index: index summary
            % - Parameter options.normalization: active normalization; omitted uses the EVP default
            % - Parameter options.metadata: additional metadata
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.N2Function: buoyancy frequency squared function
            % - Returns basisSet: initialized basis set
            arguments
                options.solver = []
                options.evp = []
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

            self.solver = options.solver;
            self.evp = options.evp;
            self.nativeModes = options.nativeModes;
            self.eigenvalues = reshape(options.eigenvalues,1,[]);
            self.h = reshape(options.h,1,[]);
            nModes = max(size(options.nativeModes,2), length(self.h));
            self.modeNumber = IMBasisSet.resolveModeNumber(options.modeNumber, nModes);
            self.index = options.index;
            self.normalization = IMBasisSet.resolveDefaultNormalization(options.normalization, options.evp);
            self.metadata = options.metadata;
            self.zDomain = options.zDomain;
            self.N2Function = options.N2Function;

            if ~isempty(self.solver)
                if any(isnan(self.zDomain))
                    self.zDomain = self.solver.zDomain;
                end
                if isempty(self.N2Function)
                    self.N2Function = @(z) self.solver.N2(z);
                end
            end
        end

        function G = G(self, z, options)
            % Evaluate $$G_j(z)$$ modes on a physical grid.
            %
            % If the EVP formulation solves `G`, this returns the solved
            % modes. If the EVP formulation solves `F`, this evaluates
            % $$G_j=-gN^{-2}\partial_zF_j$$.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: G = G(basisSet,z,options)
            % - Parameter z: physical evaluation points
            % - Parameter options.normalization: normalization to apply
            % - Returns G: evaluated `G` modes
            arguments
                self IMBasisSet
                z (:,1) double
                options.normalization = self.normalization
            end

            rawValues = self.rawVariable("G", z);
            factors = self.normalizationFactors(options.normalization);
            G = rawValues ./ factors;
        end

        function F = F(self, z, options)
            % Evaluate $$F_j(z)$$ modes on a physical grid.
            %
            % If the EVP formulation solves `F`, this returns the solved
            % modes. If the EVP formulation solves `G`, this evaluates
            % $$F_j=h_j\partial_zG_j$$.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: F = F(basisSet,z,options)
            % - Parameter z: physical evaluation points
            % - Parameter options.normalization: normalization to apply
            % - Returns F: evaluated `F` modes
            arguments
                self IMBasisSet
                z (:,1) double
                options.normalization = self.normalization
            end

            rawValues = self.rawVariable("F", z);
            factors = self.normalizationFactors(options.normalization);
            F = rawValues ./ factors;
        end

        function values = evaluate(self, variable, z, options)
            % Evaluate `F` or `G` on a physical grid.
            %
            % This method is a validated dispatcher around the canonical
            % `F(z)` and `G(z)` methods.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = evaluate(basisSet,variable,z,options)
            % - Parameter variable: `"F"` or `"G"`
            % - Parameter z: physical evaluation points
            % - Parameter options.normalization: normalization to apply
            % - Returns values: evaluated variable matrix
            arguments
                self IMBasisSet
                variable {mustBeTextScalar}
                z (:,1) double
                options.normalization = self.normalization
            end

            switch IMBasisSet.validateVariable(variable)
                case "G"
                    values = self.G(z, normalization=options.normalization);
                case "F"
                    values = self.F(z, normalization=options.normalization);
            end
        end

        function values = evaluateAll(self, z)
            % Evaluate `G` and `F` on a physical grid.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = evaluateAll(basisSet,z)
            % - Parameter z: physical evaluation points
            % - Returns values: structure with fields `G` and `F`
            values.G = self.G(z);
            values.F = self.F(z);
        end

        function factors = normalizationFactors(self, normalization)
            % Return factors for the requested normalization.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: factors = normalizationFactors(basisSet,normalization)
            % - Parameter normalization: normalization convention
            % - Returns factors: row vector of normalization factors
            name = char(self.normalizationName(normalization));
            if isempty(self.evp) || ~isfield(self.evp.normalizations, name)
                error("IMBasisSet:UnsupportedNormalization", ...
                    "The EVP does not define a ""%s"" normalization.", name);
            end
            normalizeMode = self.evp.normalizations.(name);
            nModes = max(length(self.h), size(self.nativeModes,2));
            factors = zeros(1,nModes);
            for iMode = 1:nModes
                factors(iMode) = normalizeMode(self, iMode);
            end
            factors = abs(factors);
            factors(factors == 0 | ~isfinite(factors)) = 1;
        end

        function modes = normalizedNativeModes(self, normalization)
            % Return native modes scaled by a normalization convention.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: modes = normalizedNativeModes(basisSet,normalization)
            % - Parameter normalization: normalization convention
            % - Returns modes: scaled native modes
            factors = self.normalizationFactors(normalization);
            modes = self.nativeModes ./ factors;
        end

        function gram = gramMatrix(self, variable)
            % Return the full-depth Gram matrix for a variable.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = gramMatrix(basisSet,variable)
            % - Parameter variable: variable name
            % - Returns gram: full-depth Gram matrix
            gram = self.partialGramMatrix(variable, self.zDomain(1), self.zDomain(2));
        end

        function gram = partialGramMatrix(self, variable, zMin, zMax)
            % Return a partial-depth Gram matrix.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = partialGramMatrix(basisSet,variable,zMin,zMax)
            % - Parameter variable: variable name
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns gram: partial-depth Gram matrix
            gram = self.variableGramMatrix(variable, [zMin zMax]);
        end

        function windowModes = partialWindowModes(self, variable, zMin, zMax)
            % Diagonalize partial-depth energy in existing coefficient coordinates.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: windowModes = partialWindowModes(basisSet,variable,zMin,zMax)
            % - Parameter variable: variable name
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns windowModes: window-mode decomposition structure
            gram = self.partialGramMatrix(variable, zMin, zMax);
            gram = 0.5*(gram + gram.');
            [R, D] = eig(gram);
            [eigenvalues, sortIndex] = sort(diag(D), "descend");
            windowModes.rotation = R(:,sortIndex);
            windowModes.eigenvalues = eigenvalues(:).';
            windowModes.gramMatrix = gram;
        end

        function transform = nativeTransform(self, varargin)
            % Build a native transform.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: transform = nativeTransform(basisSet,varargin)
            % - Returns transform: native transform
            self.unsupported("nativeTransform");
        end

        function projection = observationProjection(self, varargin)
            % Build an observation projection.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: projection = observationProjection(basisSet,varargin)
            % - Returns projection: observation projection
            self.unsupported("observationProjection");
        end

        function spectrum = spectrum(self, coefficients, options)
            % Compute a modal spectrum from coefficients.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = spectrum(basisSet,coefficients,options)
            % - Parameter coefficients: modal coefficients
            % - Parameter options.variable: variable name
            % - Returns spectrum: modal spectrum
            arguments
                self IMBasisSet
                coefficients (:,1) double
                options.variable {mustBeTextScalar} = self.evp.formulation
            end

            spectrum = self.crossSpectrum(coefficients, coefficients, variable=options.variable);
        end

        function spectrum = crossSpectrum(self, coefficientsA, coefficientsB, options)
            % Compute a modal cross-spectrum from coefficients.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB,options)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Parameter options.variable: variable name
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMBasisSet
                coefficientsA (:,1) double
                coefficientsB (:,1) double
                options.variable {mustBeTextScalar} = self.evp.formulation
            end

            gram = self.gramMatrix(options.variable);
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
        end

        function values = N2(self, z)
            % Evaluate buoyancy frequency squared for this basis set.
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
            % Evaluate $$\partial_z\log N^2$$ for this basis set.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = dzLogN2(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: vertical derivative of `log(N2)`
            if ~isempty(self.solver)
                values = self.solver.dzLogN2(z);
                return;
            end
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
            % Orient modes with a deterministic physical sign convention.
            %
            % Modes are signed so that the surface `F` trace is positive
            % whenever it is nonzero. If the surface `F` trace is zero, the
            % convention falls back to `G`, first at the surface and then at
            % its largest resolved amplitude.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: basisSet = orientModeSigns(basisSet)
            % - Returns basisSet: basis set with oriented native modes
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
                [reference, tolerance] = self.surfaceSignReference(FSurface, FGrid, iMode);
                if ~isfinite(reference) || abs(reference) <= tolerance
                    [reference, tolerance] = self.surfaceOrDominantSignReference(GSurface, GGrid, iMode);
                end
                if ~isfinite(reference) || abs(reference) <= tolerance
                    reference = self.dominantNativeValue(iMode);
                end
                if reference < 0
                    signs(iMode) = -1;
                end
            end
            self.nativeModes = self.nativeModes .* signs;
        end

        function factor = innerProductNormFactor(self, variable, iMode)
            % Return the norm factor from a variable's declared inner product.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = innerProductNormFactor(basisSet,variable,iMode)
            % - Parameter variable: variable name
            % - Parameter iMode: mode index
            % - Returns factor: divisor that makes the variable unit norm
            innerWeight = self.variableInteriorWeight(variable);
            innerProductTerms = self.variableInnerProductTerms(variable);
            [surfaceWeight, bottomWeight, remainingTerms] = self.scalarEndpointWeights(variable, innerProductTerms);
            if isempty(remainingTerms)
                factor = self.weightedNormFactor(variable, iMode, innerWeight, surfaceWeight, bottomWeight);
            else
                factor = self.weightedNormFactorWithInnerProductTerms(variable, iMode, innerWeight, innerProductTerms);
            end
        end

        function factor = weightedNormFactor(self, variable, iMode, innerWeight, surfaceWeight, bottomWeight)
            % Return the norm factor from an explicit quadratic form.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = weightedNormFactor(basisSet,variable,iMode,innerWeight,surfaceWeight,bottomWeight)
            % - Parameter variable: variable name
            % - Parameter iMode: mode index
            % - Parameter innerWeight: interior weight coefficient
            % - Parameter surfaceWeight: surface endpoint weight
            % - Parameter bottomWeight: bottom endpoint weight
            % - Returns factor: divisor from the requested quadratic form
            z = self.innerProductGrid(self.zDomain);
            values = self.rawVariable(variable, z);
            weight = IMOperator.evaluateCoefficient(innerWeight, z, self.context());
            integrand = weight(:).*values(:,iMode).*values(:,iMode);
            normValue = self.integrateInnerProduct(z, integrand, self.zDomain);
            innerProductTerms = IMBasisSet.scalarInnerProductTerms(variable, surfaceWeight, bottomWeight);
            normValue = normValue + self.endpointQuadraticValue(innerProductTerms, iMode, iMode, false);
            factor = sqrt(abs(normValue));
        end

        function factor = weightedNormFactorWithBoundaryTerms(self, variable, iMode, innerWeight)
            % Return a norm factor using inner-product terms from resolved boundaries.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = weightedNormFactorWithBoundaryTerms(basisSet,variable,iMode,innerWeight)
            % - Parameter variable: variable name
            % - Parameter iMode: mode index
            % - Parameter innerWeight: interior weight coefficient
            % - Returns factor: divisor from the boundary-aware quadratic form
            innerProductTerms = self.variableInnerProductTerms(variable);
            factor = self.weightedNormFactorWithInnerProductTerms(variable, iMode, innerWeight, innerProductTerms);
        end

        function factor = weightedInteriorNormFactor(self, variable, iMode, innerWeight)
            % Return a norm factor from only an interior weighted integral.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = weightedInteriorNormFactor(basisSet,variable,iMode,innerWeight)
            % - Parameter variable: variable name
            % - Parameter iMode: mode index
            % - Parameter innerWeight: interior weight coefficient
            % - Returns factor: divisor from the interior quadratic form
            factor = self.weightedNormFactorWithInnerProductTerms( ...
                variable, iMode, innerWeight, IMBasisSet.emptyInnerProductTerms());
        end

        function factor = weightedNormFactorWithInnerProductTerms(self, variable, iMode, innerWeight, innerProductTerms)
            % Return a norm factor from boundary inner-product terms.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = weightedNormFactorWithInnerProductTerms(basisSet,variable,iMode,innerWeight,innerProductTerms)
            % - Parameter variable: variable name
            % - Parameter iMode: mode index
            % - Parameter innerWeight: interior weight coefficient
            % - Parameter innerProductTerms: boundary inner-product terms
            % - Returns factor: divisor from the requested quadratic form
            z = self.innerProductGrid(self.zDomain);
            values = self.rawVariable(variable, z);
            weight = IMOperator.evaluateCoefficient(innerWeight, z, self.context());
            integrand = weight(:).*values(:,iMode).*values(:,iMode);
            normValue = self.integrateInnerProduct(z, integrand, self.zDomain);
            normValue = normValue + self.endpointQuadraticValue(innerProductTerms, iMode, iMode, false);
            factor = sqrt(abs(normValue));
        end

        function factor = geostrophicNormFactor(self, iMode)
            % Return the hydrostatic geostrophic normalization factor.
            %
            % Baroclinic modes are normalized by the `G` metric
            % $$g^{-1}\int N^2G^2\,dz$$. The `F`-solved null mode has
            % $$F_0(z)=1$$ and $$G_0(z)=0$$, so it is normalized by its RMS
            % `F` amplitude and retains Parseval scale $$D$$.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = geostrophicNormFactor(basisSet,iMode)
            % - Parameter iMode: mode index
            % - Returns factor: divisor for geostrophic normalization
            modeNumber = self.modeNumber(iMode);
            if modeNumber < 0 && self.hasUnknownBoundaryInnerProduct()
                error("IMBasisSet:UnsupportedNormalization", ...
                    "Geostrophic normalization for boundary mode %d requires known boundary inner-product metadata.", modeNumber);
            end

            fFactor = self.weightedInteriorNormFactor("F", iMode, @(z,ctx) ones(size(z)));
            if modeNumber == 0
                gFactor = self.weightedNormFactorWithBoundaryTerms("G", iMode, @(z,ctx) ctx.N2(z)/ctx.g);
                self.validateHydrostaticFNullMode(iMode, gFactor, fFactor);
                factor = fFactor/sqrt(diff(self.zDomain));
                return;
            end

            factor = self.weightedNormFactorWithBoundaryTerms("G", iMode, @(z,ctx) ctx.N2(z)/ctx.g);
        end

        function factor = maxAbsFactor(self, variable, iMode)
            % Return the maximum-amplitude normalization factor.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = maxAbsFactor(basisSet,variable,iMode)
            % - Parameter variable: variable name
            % - Parameter iMode: mode index
            % - Returns factor: maximum absolute variable amplitude
            z = self.integrationGrid(self.zDomain);
            values = self.rawVariable(variable, z);
            factor = max(abs(values(:,iMode)));
        end

        function factor = surfacePressureNormFactor(self, iMode)
            % Return the surface-pressure normalization factor.
            %
            % Surface-pressure normalization divides by the raw surface
            % trace of `F`,
            % $$A_j=F_j^\mathrm{raw}(z_\mathrm{surface})$$, so that
            % $$F_j^\mathrm{surfacePressure}(z)=F_j^\mathrm{raw}(z)/A_j$$
            % and $$F_j^\mathrm{surfacePressure}(z_\mathrm{surface})=1$$.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = surfacePressureNormFactor(basisSet,iMode)
            % - Parameter iMode: mode index
            % - Returns factor: raw surface `F` trace
            zSurface = self.zDomain(2);
            surfaceValues = self.rawVariable("F", zSurface);
            surfaceValue = surfaceValues(1,iMode);
            z = self.integrationGrid(self.zDomain);
            values = self.rawVariable("F", z);
            scale = max([1 abs(surfaceValue) max(abs(values(:,iMode)))]);
            tolerance = 1e-10*scale;
            if ~isfinite(surfaceValue) || abs(surfaceValue) <= tolerance
                error("IMBasisSet:UnsupportedNormalization", ...
                    "Surface-pressure normalization requires a finite nonzero surface F trace for mode %d.", iMode);
            end
            factor = surfaceValue;
        end
    end

    methods (Access = protected)
        function values = rawVariable(self, variable, z)
            variable = IMBasisSet.validateVariable(variable);
            if variable == self.evp.formulation
                self.requireSolver("evaluate solved variables");
                values = self.solver.evaluateNativeModes(self.nativeModes, z);
                return;
            end

            self.requireSolver("evaluate diagnostic variables");
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

        function z = integrationGrid(self, zBounds)
            if ~isempty(self.solver)
                nGrid = max(256, 4*self.solver.nEVP);
            else
                nGrid = max(256, 4*length(self.h));
            end
            z = linspace(min(zBounds), max(zBounds), nGrid).';
        end

        function gram = variableGramMatrix(self, variable, zBounds)
            z = self.innerProductGrid(zBounds);
            values = self.evaluate(variable, z);
            innerWeight = self.variableInteriorWeight(variable);
            innerProductTerms = self.variableInnerProductTerms(variable);
            weight = IMOperator.evaluateCoefficient(innerWeight, z, self.context());
            gram = zeros(size(values,2), size(values,2));
            for iMode = 1:size(values,2)
                for jMode = iMode:size(values,2)
                    integrand = weight(:).*values(:,iMode).*values(:,jMode);
                    value = self.integrateInnerProduct(z, integrand, zBounds);
                    value = value + self.endpointQuadraticValue(innerProductTerms, iMode, jMode, true, zBounds);
                    gram(iMode,jMode) = value;
                    gram(jMode,iMode) = value;
                end
            end
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

        function innerWeight = variableInteriorWeight(self, variable)
            variable = IMBasisSet.validateVariable(variable);
            fieldName = char(variable);
            if ~isempty(self.evp) && isfield(self.evp.innerWeights, fieldName)
                innerWeight = self.evp.innerWeights.(fieldName);
            elseif variable == "G"
                innerWeight = @(z,ctx) ctx.N2(z)/ctx.g;
            else
                innerWeight = @(z,ctx) ones(size(z));
            end
        end

        function innerProductTerms = variableInnerProductTerms(self, variable)
            innerProductTerms = IMBasisSet.emptyInnerProductTerms();
            if isempty(self.evp) || ~isprop(self.evp, "boundaryConditions")
                return;
            end

            variable = IMBasisSet.validateVariable(variable);
            for iBoundary = 1:length(self.evp.boundaryConditions)
                terms = self.evp.boundaryConditions(iBoundary).innerProductTerms;
                for iTerm = 1:length(terms)
                    if string(terms(iTerm).innerProductVariable) == variable
                        innerProductTerms(end+1,1) = terms(iTerm);
                    end
                end
            end
        end

        function value = endpointQuadraticValue(self, innerProductTerms, iMode, jMode, useNormalized, zBounds)
            value = 0;
            if nargin < 6
                zBounds = self.zDomain;
            end
            for iTerm = 1:length(innerProductTerms)
                value = value + self.innerProductTermValue(innerProductTerms(iTerm), iMode, jMode, useNormalized, zBounds);
            end
        end

        function value = innerProductTermValue(self, term, iMode, jMode, useNormalized, zBounds)
            value = 0;
            switch string(term.location)
                case "surface"
                    zEndpoint = self.zDomain(2);
                case "bottom"
                    zEndpoint = self.zDomain(1);
                otherwise
                    error("IMBasisSet:InvalidEndpointLocation", ...
                        "Unknown endpoint location ""%s"".", string(term.location));
            end
            if ~self.boundsIncludeEndpoint(zBounds, zEndpoint)
                return;
            end

            coefficient = self.endpointCoefficient(term.coefficient, zEndpoint);
            leftValues = self.traceValues(term.leftTrace, zEndpoint, useNormalized);
            rightValues = self.traceValues(term.rightTrace, zEndpoint, useNormalized);
            value = coefficient*leftValues(:,iMode).*rightValues(:,jMode);
        end

        function values = traceValues(self, trace, z, useNormalized)
            variable = IMBasisSet.validateVariable(trace.variable);
            derivativeOrder = trace.derivativeOrder;
            if derivativeOrder == 0
                if useNormalized
                    values = self.evaluate(variable, z);
                else
                    values = self.rawVariable(variable, z);
                end
                return;
            end
            if derivativeOrder ~= 1
                error("IMBasisSet:UnsupportedTrace", ...
                    "Endpoint traces currently support only variable values and first derivatives.");
            end
            if isempty(self.evp) || variable ~= self.evp.formulation
                self.unsupported("evaluate derivative traces for " + variable);
            end
            self.requireSolver("evaluate derivative traces");
            values = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, derivativeOrder);
            if useNormalized
                values = values ./ self.normalizationFactors(self.normalization);
            end
        end

        function coefficient = endpointCoefficient(self, coefficient, z)
            if isa(coefficient, "function_handle")
                context = self.context();
                try
                    coefficient = coefficient(context);
                catch
                    try
                        coefficient = coefficient(z, context);
                    catch
                        coefficient = coefficient(z);
                    end
                end
            end
            if ~isscalar(coefficient)
                coefficient = coefficient(1);
            end
        end

        function [surfaceWeight, bottomWeight, remainingTerms] = scalarEndpointWeights(~, variable, innerProductTerms)
            variable = IMBasisSet.validateVariable(variable);
            surfaceWeight = 0;
            bottomWeight = 0;
            remainingTerms = IMBasisSet.emptyInnerProductTerms();
            for iTerm = 1:length(innerProductTerms)
                term = innerProductTerms(iTerm);
                leftTrace = term.leftTrace;
                rightTrace = term.rightTrace;
                isScalarTerm = isnumeric(term.coefficient) && isscalar(term.coefficient) ...
                    && string(leftTrace.variable) == variable && string(rightTrace.variable) == variable ...
                    && leftTrace.derivativeOrder == 0 && rightTrace.derivativeOrder == 0;
                if isScalarTerm && string(term.location) == "surface"
                    surfaceWeight = surfaceWeight + term.coefficient;
                elseif isScalarTerm && string(term.location) == "bottom"
                    bottomWeight = bottomWeight + term.coefficient;
                else
                    remainingTerms(end+1,1) = term;
                end
            end
        end

        function tf = hasUnknownBoundaryInnerProduct(self)
            tf = false;
            if isempty(self.evp) || ~isprop(self.evp, "boundaryConditions")
                return;
            end
            if isempty(self.evp.boundaryConditions)
                return;
            end
            tf = any(~[self.evp.boundaryConditions.hasKnownInnerProductTerms]);
        end

        function validateHydrostaticFNullMode(self, iMode, gFactor, fFactor)
            if isempty(self.evp) || self.evp.formulation ~= "F"
                error("IMBasisSet:UnsupportedNormalization", ...
                    "Mode number 0 geostrophic normalization is only defined for hydrostatic F modes.");
            end
            if self.evp.name ~= "hydrostaticFModes"
                error("IMBasisSet:UnsupportedNormalization", ...
                    "Mode number 0 geostrophic normalization is only defined for hydrostatic F modes.");
            end
            if iMode > length(self.eigenvalues)
                return;
            end

            positiveEigenvalues = self.eigenvalues(self.eigenvalues > 0 & isfinite(self.eigenvalues));
            if isempty(positiveEigenvalues)
                eigenvalueTolerance = 1e-10;
            else
                eigenvalueTolerance = 1e-8*max(1,min(positiveEigenvalues));
            end
            gTolerance = 1e-8*max(1,fFactor);
            if abs(self.eigenvalues(iMode)) > eigenvalueTolerance || gFactor > gTolerance
                warning("IMBasisSet:ModeNumberValidation", ...
                    "Mode number 0 is expected to be a hydrostatic F null mode with lambda near zero and G near zero.");
            end
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

        function [reference, tolerance] = surfaceSignReference(~, surfaceValues, gridValues, iMode)
            reference = NaN;
            scale = 0;
            if ~isempty(surfaceValues)
                reference = surfaceValues(1,iMode);
                scale = max(scale, abs(reference));
            end
            if ~isempty(gridValues)
                scale = max(scale, max(abs(gridValues(:,iMode))));
            end
            tolerance = 1e-10*max(1,scale);
        end

        function [reference, tolerance] = surfaceOrDominantSignReference(self, surfaceValues, gridValues, iMode)
            [reference, tolerance] = self.surfaceSignReference(surfaceValues, gridValues, iMode);
            if isfinite(reference) && abs(reference) > tolerance
                return;
            end
            if isempty(gridValues)
                return;
            end

            [scale, index] = max(abs(gridValues(:,iMode)));
            reference = gridValues(index,iMode);
            tolerance = 1e-10*max(1,scale);
        end

        function value = dominantNativeValue(self, iMode)
            [~, index] = max(abs(self.nativeModes(:,iMode)));
            value = self.nativeModes(index,iMode);
        end

        function context = context(self)
            context.N2 = @(z) self.N2(z);
            context.dzLogN2 = @(z) self.dzLogN2(z);
            context.zDomain = self.zDomain;
            context.coordinateKind = "basisSet";
            if ~isempty(self.evp)
                context.f0 = self.evp.f0;
                context.g = self.evp.g;
            else
                context.f0 = 0;
                context.g = 9.81;
            end
            if ~isempty(self.solver)
                solverContext = self.solver.context();
                if isfield(solverContext, "coordinateKind")
                    context.coordinateKind = solverContext.coordinateKind;
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
            % Resolve an omitted normalization from the EVP default.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: normalization = IMBasisSet.resolveDefaultNormalization(requestedNormalization,evp)
            % - Parameter requestedNormalization: caller-specified normalization or empty
            % - Parameter evp: EVP descriptor
            % - Returns normalization: resolved normalization
            if ~isempty(requestedNormalization)
                normalization = requestedNormalization;
                return;
            end

            if ~isempty(evp) && isprop(evp, "defaultNormalization") && ~isempty(evp.defaultNormalization)
                normalization = evp.defaultNormalization;
                return;
            end

            normalization = Normalization.kConstant;
        end

        function variable = validateVariable(variable)
            variable = string(variable);
            if variable ~= "F" && variable ~= "G"
                error("IMBasisSet:InvalidVariable", ...
                    "variable must be ""F"" or ""G"".");
            end
        end
    end

    methods (Static)
        function terms = emptyInnerProductTerms()
            % Create an empty boundary inner-product-term structure.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: terms = IMBasisSet.emptyInnerProductTerms()
            % - Returns terms: empty boundary inner-product-term structure
            terms = struct("innerProductVariable", {}, "location", {}, "coefficient", {}, ...
                "leftTrace", {}, "rightTrace", {});
        end

        function terms = scalarInnerProductTerms(variable, surfaceWeight, bottomWeight)
            % Convert scalar endpoint weights to boundary inner-product terms.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: terms = IMBasisSet.scalarInnerProductTerms(variable,surfaceWeight,bottomWeight)
            % - Parameter variable: variable name
            % - Parameter surfaceWeight: surface endpoint weight
            % - Parameter bottomWeight: bottom endpoint weight
            % - Returns terms: boundary inner-product terms
            variable = IMBasisSet.validateVariable(variable);
            terms = IMBasisSet.emptyInnerProductTerms();
            if surfaceWeight ~= 0
                terms(end+1,1) = IMBoundary.innerProductTerm(variable, "surface", surfaceWeight, ...
                    IMBoundary.trace(variable), IMBoundary.trace(variable));
            end
            if bottomWeight ~= 0
                terms(end+1,1) = IMBoundary.innerProductTerm(variable, "bottom", bottomWeight, ...
                    IMBoundary.trace(variable), IMBoundary.trace(variable));
            end
        end

        function basisSet = constantStratification(options)
            % Create an analytical constant-stratification basis set.
            %
            % The returned basis set evaluates exact constant-stratification
            % `F` and `G` modes without solving a numerical EVP.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet.constantStratification(options)
            % - Parameter options.evp: supported wave-mode or hydrostatic eigenvalue-problem descriptor
            % - Parameter options.N0: constant buoyancy frequency
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.nModes: number of modes
            % - Parameter options.normalization: active normalization; omitted uses the EVP default
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: analytical constant-stratification basis set
            arguments
                options.evp IMEigenvalueProblem = IMEigenvalueProblem.hydrostaticGModes()
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
            % The returned basis set evaluates exact rigid-endpoint
            % `G`-formulation modes for
            % $$N^2(z)=N_0^2e^{2z/b}$$ without solving a numerical EVP.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet.exponentialStratification(options)
            % - Parameter options.evp: supported rigid-endpoint `G` eigenvalue-problem descriptor
            % - Parameter options.N0: surface buoyancy frequency
            % - Parameter options.b: exponential e-folding depth
            % - Parameter options.zDomain: physical vertical domain with surface at zero
            % - Parameter options.nModes: number of modes
            % - Parameter options.normalization: active normalization; omitted uses the EVP default
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: analytical exponential-stratification basis set
            arguments
                options.evp IMEigenvalueProblem = IMEigenvalueProblem.hydrostaticGModes()
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

        function basisSet = wkbApproximation(options)
            % Report that WKB analytical bases are deferred.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet.wkbApproximation(options)
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: not returned because this factory is not implemented
            arguments
                options.metadata struct = struct()
            end

            error("IMBasisSet:AnalyticalBasisNotImplemented", ...
                "The WKB analytical basis set is not implemented yet.");
        end
    end
end
