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
        % new EVP.
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

        % Mode index coordinate.
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

        % Coriolis parameter.
        %
        % - Topic: Inspect basis sets
        f0

        % Gravitational acceleration.
        %
        % - Topic: Inspect basis sets
        g
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
            % - Parameter options.index: index summary
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.N2Function: buoyancy frequency squared function
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Returns basisSet: initialized basis set
            arguments
                options.solver = []
                options.evp = []
                options.nativeModes double = zeros(0,0)
                options.eigenvalues double = zeros(1,0)
                options.h double = zeros(1,0)
                options.index struct = struct()
                options.normalization = Normalization.kConstant
                options.metadata struct = struct()
                options.zDomain (1,2) double = [NaN NaN]
                options.N2Function = []
                options.f0 (1,1) double = NaN
                options.g (1,1) double = NaN
            end

            self.solver = options.solver;
            self.evp = options.evp;
            self.nativeModes = options.nativeModes;
            self.eigenvalues = reshape(options.eigenvalues,1,[]);
            self.h = reshape(options.h,1,[]);
            self.modeIndex = 1:max(size(options.nativeModes,2), length(self.h));
            self.index = options.index;
            self.normalization = options.normalization;
            self.metadata = options.metadata;
            self.zDomain = options.zDomain;
            self.N2Function = options.N2Function;
            self.f0 = options.f0;
            self.g = options.g;

            if ~isempty(self.solver)
                if any(isnan(self.zDomain))
                    self.zDomain = self.solver.zDomain;
                end
                if isempty(self.N2Function)
                    self.N2Function = @(z) self.solver.N2(z);
                end
                if isnan(self.f0)
                    self.f0 = self.solver.f0;
                end
                if isnan(self.g)
                    self.g = self.solver.g;
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
    end

    methods (Hidden)
        function factor = innerProductNormFactor(self, component, iMode)
            % Return the norm factor from a component's declared inner product.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: factor = innerProductNormFactor(basisSet,component,iMode)
            % - Parameter component: component name
            % - Parameter iMode: mode index
            % - Returns factor: divisor that makes the component unit norm
            [innerWeight, surfaceWeight, bottomWeight] = self.componentInnerProductWeights(component);
            factor = self.weightedNormFactor(component, iMode, innerWeight, surfaceWeight, bottomWeight);
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
            normValue = normValue + self.endpointQuadraticValue(component, iMode, iMode, surfaceWeight, bottomWeight);
            factor = sqrt(abs(normValue));
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
            values = metadata.operator.evaluate(self.solver, self.nativeModes, z);
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
            [innerWeight, surfaceWeight, bottomWeight] = self.componentInnerProductWeights(component);
            weight = IMOperator.evaluateCoefficient(innerWeight, z, self.context());
            gram = zeros(size(values,2), size(values,2));
            for iMode = 1:size(values,2)
                for jMode = iMode:size(values,2)
                    integrand = weight(:).*values(:,iMode).*values(:,jMode);
                    value = self.integrateInnerProduct(z, integrand, zBounds);
                    value = value + self.endpointGramValue(component, iMode, jMode, ...
                        self.activeSurfaceWeight(zBounds, surfaceWeight), ...
                        self.activeBottomWeight(zBounds, bottomWeight));
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

        function [innerWeight, surfaceWeight, bottomWeight] = componentInnerProductWeights(self, component)
            metadata = self.componentMetadata(component);
            if isfield(metadata, "innerWeight")
                innerWeight = metadata.innerWeight;
            elseif string(component) == "G"
                innerWeight = @(z,ctx) ctx.N2(z)/ctx.g;
            else
                innerWeight = @(z,ctx) ones(size(z));
            end
            surfaceWeight = 0;
            bottomWeight = 0;
            if isfield(metadata, "surfaceWeight")
                surfaceWeight = metadata.surfaceWeight;
            end
            if isfield(metadata, "bottomWeight")
                bottomWeight = metadata.bottomWeight;
            end
        end

        function value = endpointQuadraticValue(self, component, iMode, jMode, surfaceWeight, bottomWeight)
            value = 0;
            if surfaceWeight ~= 0
                surfaceValues = self.rawComponent(component, self.zDomain(2));
                value = value + surfaceWeight*surfaceValues(:,iMode).*surfaceValues(:,jMode);
            end
            if bottomWeight ~= 0
                bottomValues = self.rawComponent(component, self.zDomain(1));
                value = value + bottomWeight*bottomValues(:,iMode).*bottomValues(:,jMode);
            end
        end

        function value = endpointGramValue(self, component, iMode, jMode, surfaceWeight, bottomWeight)
            value = 0;
            if surfaceWeight ~= 0
                surfaceValues = self.evaluate(component, self.zDomain(2));
                value = value + surfaceWeight*surfaceValues(:,iMode).*surfaceValues(:,jMode);
            end
            if bottomWeight ~= 0
                bottomValues = self.evaluate(component, self.zDomain(1));
                value = value + bottomWeight*bottomValues(:,iMode).*bottomValues(:,jMode);
            end
        end

        function weight = activeSurfaceWeight(self, zBounds, surfaceWeight)
            weight = 0;
            if self.boundsIncludeEndpoint(zBounds, self.zDomain(2))
                weight = surfaceWeight;
            end
        end

        function weight = activeBottomWeight(self, zBounds, bottomWeight)
            weight = 0;
            if self.boundsIncludeEndpoint(zBounds, self.zDomain(1))
                weight = bottomWeight;
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

        function context = context(self)
            context.N2 = @(z) self.N2(z);
            context.f0 = self.f0;
            context.g = self.g;
            context.zDomain = self.zDomain;
            context.parameters = struct();
            if ~isempty(self.evp)
                context.parameters = self.evp.parameters;
                names = fieldnames(self.evp.parameters);
                for iName = 1:length(names)
                    context.(names{iName}) = self.evp.parameters.(names{iName});
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

    methods (Static)
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
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: analytical constant-stratification basis set
            arguments
                options.evp IMEigenvalueProblem = IMEigenvalueProblem.hydrostaticGModes()
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.zDomain (1,2) double = [-1 0]
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.f0 (1,1) double = NaN
                options.g (1,1) double = NaN
                options.normalization = Normalization.kConstant
                options.metadata struct = struct()
            end

            basisSet = IMBasisSetConstantStratification(evp=options.evp, N0=options.N0, ...
                zDomain=options.zDomain, nModes=options.nModes, f0=options.f0, g=options.g, ...
                normalization=options.normalization, metadata=options.metadata);
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
