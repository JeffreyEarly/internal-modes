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
    % G = basisSet.evaluate("G", linspace(-1000,0,64).');
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

        % Scientific mode labels.
        %
        % Negative labels identify boundary modes, zero labels identify
        % barotropic or null modes, and positive labels identify baroclinic
        % modes.
        %
        % - Topic: Inspect basis sets
        modeIndex

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
            % - Parameter options.modeIndex: scientific mode labels
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
                options.modeIndex double = []
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
            self.modeIndex = IMBasisSet.resolveModeIndex(options.modeIndex, nModes);
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

        function values = evaluate(self, component, z, options)
            % Evaluate a component on a physical grid.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = evaluate(basisSet,component,z,options)
            % - Parameter component: component name
            % - Parameter z: physical evaluation points
            % - Parameter options.normalization: normalization to apply
            % - Returns values: evaluated component matrix
            arguments
                self IMBasisSet
                component {mustBeTextScalar}
                z (:,1) double
                options.normalization = self.normalization
            end

            rawValues = self.rawComponent(component, z);
            factors = self.normalizationFactors(options.normalization);
            values = rawValues ./ factors;
        end

        function values = evaluateAll(self, z)
            % Evaluate all EVP-declared components on a physical grid.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = evaluateAll(basisSet,z)
            % - Parameter z: physical evaluation points
            % - Returns values: structure of evaluated components
            values = struct();
            componentNames = fieldnames(self.evp.components);
            for iComponent = 1:length(componentNames)
                values.(componentNames{iComponent}) = self.evaluate(componentNames{iComponent}, z);
            end
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

        function gram = gramMatrix(self, component)
            % Return the full-depth Gram matrix for a component.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = gramMatrix(basisSet,component)
            % - Parameter component: component name
            % - Returns gram: full-depth Gram matrix
            gram = self.partialGramMatrix(component, self.zDomain(1), self.zDomain(2));
        end

        function gram = partialGramMatrix(self, component, zMin, zMax)
            % Return a partial-depth Gram matrix.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = partialGramMatrix(basisSet,component,zMin,zMax)
            % - Parameter component: component name
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns gram: partial-depth Gram matrix
            gram = self.componentGramMatrix(component, [zMin zMax]);
        end

        function windowModes = partialWindowModes(self, component, zMin, zMax)
            % Diagonalize partial-depth energy in existing coefficient coordinates.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: windowModes = partialWindowModes(basisSet,component,zMin,zMax)
            % - Parameter component: component name
            % - Parameter zMin: lower physical bound
            % - Parameter zMax: upper physical bound
            % - Returns windowModes: window-mode decomposition structure
            gram = self.partialGramMatrix(component, zMin, zMax);
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
            % - Parameter options.component: component name
            % - Returns spectrum: modal spectrum
            arguments
                self IMBasisSet
                coefficients (:,1) double
                options.component {mustBeTextScalar} = self.evp.primaryComponent
            end

            spectrum = self.crossSpectrum(coefficients, coefficients, component=options.component);
        end

        function spectrum = crossSpectrum(self, coefficientsA, coefficientsB, options)
            % Compute a modal cross-spectrum from coefficients.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB,options)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Parameter options.component: component name
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMBasisSet
                coefficientsA (:,1) double
                coefficientsB (:,1) double
                options.component {mustBeTextScalar} = self.evp.primaryComponent
            end

            gram = self.gramMatrix(options.component);
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
            FSurface = self.rawComponentForSign("F", zSurface);
            FGrid = self.rawComponentForSign("F", zGrid);
            GSurface = self.rawComponentForSign("G", zSurface);
            GGrid = self.rawComponentForSign("G", zGrid);
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

        function factor = innerProductNormFactor(self, component, iMode)
            % Return the norm factor from a component's declared inner product.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = innerProductNormFactor(basisSet,component,iMode)
            % - Parameter component: component name
            % - Parameter iMode: mode index
            % - Returns factor: divisor that makes the component unit norm
            innerWeight = self.componentInteriorWeight(component);
            innerProductTerms = self.componentInnerProductTerms(component);
            [surfaceWeight, bottomWeight, remainingTerms] = self.scalarEndpointWeights(component, innerProductTerms);
            if isempty(remainingTerms)
                factor = self.weightedNormFactor(component, iMode, innerWeight, surfaceWeight, bottomWeight);
            else
                factor = self.weightedNormFactorWithInnerProductTerms(component, iMode, innerWeight, innerProductTerms);
            end
        end

        function factor = weightedNormFactor(self, component, iMode, innerWeight, surfaceWeight, bottomWeight)
            % Return the norm factor from an explicit quadratic form.
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
            z = self.innerProductGrid(self.zDomain);
            values = self.rawComponent(component, z);
            weight = IMOperator.evaluateCoefficient(innerWeight, z, self.context());
            integrand = weight(:).*values(:,iMode).*values(:,iMode);
            normValue = self.integrateInnerProduct(z, integrand, self.zDomain);
            innerProductTerms = IMBasisSet.scalarInnerProductTerms(component, surfaceWeight, bottomWeight);
            normValue = normValue + self.endpointQuadraticValue(innerProductTerms, iMode, iMode, false);
            factor = sqrt(abs(normValue));
        end

        function factor = weightedNormFactorWithBoundaryTerms(self, component, iMode, innerWeight)
            % Return a norm factor using inner-product terms from resolved boundaries.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = weightedNormFactorWithBoundaryTerms(basisSet,component,iMode,innerWeight)
            % - Parameter component: component name
            % - Parameter iMode: mode index
            % - Parameter innerWeight: interior weight coefficient
            % - Returns factor: divisor from the boundary-aware quadratic form
            innerProductTerms = self.componentInnerProductTerms(component);
            factor = self.weightedNormFactorWithInnerProductTerms(component, iMode, innerWeight, innerProductTerms);
        end

        function factor = weightedInteriorNormFactor(self, component, iMode, innerWeight)
            % Return a norm factor from only an interior weighted integral.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = weightedInteriorNormFactor(basisSet,component,iMode,innerWeight)
            % - Parameter component: component name
            % - Parameter iMode: mode index
            % - Parameter innerWeight: interior weight coefficient
            % - Returns factor: divisor from the interior quadratic form
            factor = self.weightedNormFactorWithInnerProductTerms( ...
                component, iMode, innerWeight, IMBasisSet.emptyInnerProductTerms());
        end

        function factor = weightedNormFactorWithInnerProductTerms(self, component, iMode, innerWeight, innerProductTerms)
            % Return a norm factor from boundary inner-product terms.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = weightedNormFactorWithInnerProductTerms(basisSet,component,iMode,innerWeight,innerProductTerms)
            % - Parameter component: component name
            % - Parameter iMode: mode index
            % - Parameter innerWeight: interior weight coefficient
            % - Parameter innerProductTerms: boundary inner-product terms
            % - Returns factor: divisor from the requested quadratic form
            z = self.innerProductGrid(self.zDomain);
            values = self.rawComponent(component, z);
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
            % $$g^{-1}\int N^2G^2\,dz$$. The `F`-solved barotropic mode is
            % a null mode with `G=0`, so it is normalized by its RMS `F`
            % amplitude and retains Parseval scale $$D$$.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = geostrophicNormFactor(basisSet,iMode)
            % - Parameter iMode: mode index
            % - Returns factor: divisor for geostrophic normalization
            jMode = self.modeIndex(iMode);
            if jMode < 0 && self.hasUnknownBoundaryInnerProduct()
                error("IMBasisSet:UnsupportedNormalization", ...
                    "Geostrophic normalization for boundary mode %d requires known boundary inner-product metadata.", jMode);
            end

            fFactor = self.weightedInteriorNormFactor("F", iMode, @(z,ctx) ones(size(z)));
            if jMode == 0
                gFactor = self.weightedNormFactorWithBoundaryTerms("G", iMode, @(z,ctx) ctx.N2(z)/ctx.g);
                self.validateHydrostaticFBarotropicNull(iMode, gFactor, fFactor);
                factor = fFactor/sqrt(diff(self.zDomain));
                return;
            end

            factor = self.weightedNormFactorWithBoundaryTerms("G", iMode, @(z,ctx) ctx.N2(z)/ctx.g);
        end

        function factor = maxAbsFactor(self, component, iMode)
            % Return the maximum-amplitude normalization factor.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = maxAbsFactor(basisSet,component,iMode)
            % - Parameter component: component name
            % - Parameter iMode: mode index
            % - Returns factor: maximum absolute component amplitude
            z = self.integrationGrid(self.zDomain);
            values = self.rawComponent(component, z);
            factor = max(abs(values(:,iMode)));
        end
    end

    methods (Access = protected)
        function values = rawComponent(self, component, z)
            component = string(component);
            if component == self.evp.primaryComponent
                self.requireSolver("evaluate primary components");
                values = self.solver.evaluateNativeModes(self.nativeModes, z);
                return;
            end

            metadata = self.componentMetadata(component);
            if ~isfield(metadata, "role") || string(metadata.role) ~= "diagnostic"
                self.unsupported("evaluate " + component);
            end
            if ~isfield(metadata, "from") || string(metadata.from) ~= self.evp.primaryComponent
                self.unsupported("evaluate " + component + " diagnostic");
            end
            if ~isfield(metadata, "operator")
                self.unsupported("evaluate " + component + " diagnostic");
            end

            self.requireSolver("evaluate diagnostic components");
            values = metadata.operator.evaluate(self.solver, self.nativeModes, z, context=self.context());
            if isfield(metadata, "modalScale")
                values = self.applyModalScale(values, metadata.modalScale);
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

        function gram = componentGramMatrix(self, component, zBounds)
            z = self.innerProductGrid(zBounds);
            values = self.evaluate(component, z);
            innerWeight = self.componentInteriorWeight(component);
            innerProductTerms = self.componentInnerProductTerms(component);
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

        function innerWeight = componentInteriorWeight(self, component)
            metadata = self.componentMetadata(component);
            if isfield(metadata, "innerWeight")
                innerWeight = metadata.innerWeight;
            elseif string(component) == "G"
                innerWeight = @(z,ctx) ctx.N2(z)/ctx.g;
            else
                innerWeight = @(z,ctx) ones(size(z));
            end
        end

        function innerProductTerms = componentInnerProductTerms(self, component)
            innerProductTerms = IMBasisSet.emptyInnerProductTerms();
            if isempty(self.evp) || ~isprop(self.evp, "boundaryRows")
                return;
            end

            component = string(component);
            for iBoundary = 1:length(self.evp.boundaryRows)
                terms = self.evp.boundaryRows(iBoundary).innerProductTerms;
                for iTerm = 1:length(terms)
                    if string(terms(iTerm).innerProductComponent) == component
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
            component = string(trace.component);
            derivativeOrder = trace.derivativeOrder;
            if derivativeOrder == 0
                if useNormalized
                    values = self.evaluate(component, z);
                else
                    values = self.rawComponent(component, z);
                end
                return;
            end
            if derivativeOrder ~= 1
                error("IMBasisSet:UnsupportedTrace", ...
                    "Endpoint traces currently support only component values and first derivatives.");
            end
            if isempty(self.evp) || component ~= self.evp.primaryComponent
                self.unsupported("evaluate derivative traces for " + component);
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

        function [surfaceWeight, bottomWeight, remainingTerms] = scalarEndpointWeights(~, component, innerProductTerms)
            component = string(component);
            surfaceWeight = 0;
            bottomWeight = 0;
            remainingTerms = IMBasisSet.emptyInnerProductTerms();
            for iTerm = 1:length(innerProductTerms)
                term = innerProductTerms(iTerm);
                leftTrace = term.leftTrace;
                rightTrace = term.rightTrace;
                isScalarTerm = isnumeric(term.coefficient) && isscalar(term.coefficient) ...
                    && string(leftTrace.component) == component && string(rightTrace.component) == component ...
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

        function metadata = componentMetadata(self, component)
            component = char(string(component));
            if isempty(self.evp) || ~isfield(self.evp.components, component)
                self.unsupported("evaluate " + string(component));
            end
            metadata = self.evp.components.(component);
        end

        function value = hasComponent(self, component)
            value = ~isempty(self.evp) && isfield(self.evp.components, char(string(component)));
        end

        function tf = hasUnknownBoundaryInnerProduct(self)
            tf = false;
            if isempty(self.evp) || ~isprop(self.evp, "boundaryRows")
                return;
            end
            if isempty(self.evp.boundaryRows)
                return;
            end
            tf = any(~[self.evp.boundaryRows.hasKnownInnerProduct]);
        end

        function validateHydrostaticFBarotropicNull(self, iMode, gFactor, fFactor)
            if isempty(self.evp) || self.evp.primaryComponent ~= "F"
                error("IMBasisSet:UnsupportedNormalization", ...
                    "Mode index 0 geostrophic normalization is only defined for hydrostatic F modes.");
            end
            if ~isfield(self.evp.parameters, "problemType") || string(self.evp.parameters.problemType) ~= "hydrostaticFModes"
                error("IMBasisSet:UnsupportedNormalization", ...
                    "Mode index 0 geostrophic normalization is only defined for hydrostatic F modes.");
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
                warning("IMBasisSet:ModeIndexValidation", ...
                    "Mode index 0 is expected to be a hydrostatic F null mode with lambda near zero and G near zero.");
            end
        end

        function values = rawComponentForSign(self, component, z)
            values = [];
            if ~self.hasComponent(component)
                return;
            end
            try
                values = self.rawComponent(component, z);
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

        function values = applyModalScale(self, values, modalScale)
            if isstring(modalScale) || ischar(modalScale)
                switch string(modalScale)
                    case "h"
                        values = values .* self.h;
                    case "none"
                    otherwise
                        error("IMBasisSet:UnsupportedModalScale", ...
                            "Unsupported modal scale ""%s"".", string(modalScale));
                end
            elseif isnumeric(modalScale)
                values = values .* reshape(modalScale,1,[]);
            else
                error("IMBasisSet:UnsupportedModalScale", ...
                    "Unsupported modal scale type.");
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
        function modeIndex = resolveModeIndex(requestedModeIndex, nModes)
            if isempty(requestedModeIndex)
                modeIndex = 1:nModes;
                return;
            end
            modeIndex = reshape(requestedModeIndex,1,[]);
            if length(modeIndex) ~= nModes
                error("IMBasisSet:InvalidModeIndex", ...
                    "modeIndex must have one entry for each retained mode.");
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
    end

    methods (Static)
        function terms = emptyInnerProductTerms()
            % Create an empty boundary inner-product-term structure.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: terms = IMBasisSet.emptyInnerProductTerms()
            % - Returns terms: empty boundary inner-product-term structure
            terms = struct("innerProductComponent", {}, "location", {}, "coefficient", {}, ...
                "leftTrace", {}, "rightTrace", {});
        end

        function terms = scalarInnerProductTerms(component, surfaceWeight, bottomWeight)
            % Convert scalar endpoint weights to boundary inner-product terms.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: terms = IMBasisSet.scalarInnerProductTerms(component,surfaceWeight,bottomWeight)
            % - Parameter component: component name
            % - Parameter surfaceWeight: surface endpoint weight
            % - Parameter bottomWeight: bottom endpoint weight
            % - Returns terms: boundary inner-product terms
            terms = IMBasisSet.emptyInnerProductTerms();
            if surfaceWeight ~= 0
                terms(end+1,1) = IMBoundaryRow.innerProductTerm(component, "surface", surfaceWeight, ...
                    IMBoundaryRow.trace(component), IMBoundaryRow.trace(component));
            end
            if bottomWeight ~= 0
                terms(end+1,1) = IMBoundaryRow.innerProductTerm(component, "bottom", bottomWeight, ...
                    IMBoundaryRow.trace(component), IMBoundaryRow.trace(component));
            end
        end

        function basisSet = constantStratification(options)
            % Create an analytical constant-stratification basis set.
            %
            % The returned basis set evaluates exact constant-stratification
            % `G` modes without solving a numerical EVP.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet.constantStratification(options)
            % - Parameter options.evp: wave-mode or hydrostatic `G` eigenvalue-problem descriptor
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
            % Report that exponential-stratification analytical bases are deferred.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = IMBasisSet.exponentialStratification(options)
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: not returned because this factory is not implemented
            arguments
                options.metadata struct = struct()
            end

            error("IMBasisSet:AnalyticalBasisNotImplemented", ...
                "The exponential-stratification v2 analytical basis set is not implemented yet.");
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
                "The WKB v2 analytical basis set is not implemented yet.");
        end
    end
end
