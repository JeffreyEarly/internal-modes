classdef InternalModesBasisSet
    % Store native EVP solutions with basis-owned normalization.
    %
    % `InternalModesBasisSet` stores the raw native modes returned by a
    % solver. Normalization is applied lazily when modes, spectra, or Gram
    % matrices are requested.
    %
    % ```matlab
    % basisSet.normalization = Normalization.geostrophic;
    % G = basisSet.evaluate("G", linspace(-1000,0,64).');
    % ```
    %
    % - Topic: Create basis sets
    % - Topic: Evaluate basis sets
    % - Topic: Analyze Gram matrices
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesBasisSet

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
    end

    methods
        function self = InternalModesBasisSet(options)
            % Create a native solution basis set.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = InternalModesBasisSet(options)
            % - Parameter options.solver: solver reference
            % - Parameter options.evp: EVP descriptor
            % - Parameter options.nativeModes: native mode columns
            % - Parameter options.eigenvalues: retained eigenvalues
            % - Parameter options.h: equivalent depths
            % - Parameter options.index: index summary
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
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
            end

            self.solver = options.solver;
            self.evp = options.evp;
            self.nativeModes = options.nativeModes;
            self.eigenvalues = reshape(options.eigenvalues,1,[]);
            self.h = reshape(options.h,1,[]);
            self.modeIndex = 1:size(options.nativeModes,2);
            self.index = options.index;
            self.normalization = options.normalization;
            self.metadata = options.metadata;
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
                self InternalModesBasisSet
                component {mustBeTextScalar}
                z (:,1) double
                options.normalization = self.normalization
            end

            rawValues = self.rawComponent(component, z);
            factors = self.normalizationFactors(options.normalization);
            values = rawValues ./ factors;
        end

        function values = evaluateAll(self, z)
            % Evaluate all standard components on a physical grid.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: values = evaluateAll(basisSet,z)
            % - Parameter z: physical evaluation points
            % - Returns values: structure of evaluated components
            values = struct();
            values.(char(self.evp.primaryComponent)) = self.evaluate(self.evp.primaryComponent, z);
            if self.evp.primaryComponent == "G"
                values.F = self.evaluate("F", z);
            end
        end

        function factors = normalizationFactors(self, normalization)
            % Return factors for the requested normalization.
            %
            % - Topic: Evaluate basis sets
            % - Declaration: factors = normalizationFactors(basisSet,normalization)
            % - Parameter normalization: normalization convention
            % - Returns factors: row vector of normalization factors
            self.requireSolver("normalizationFactors");
            z = linspace(self.solver.zDomain(1), self.solver.zDomain(2), max(256,4*self.solver.nEVP)).';
            G = [];
            F = [];
            if self.evp.primaryComponent == "G"
                G = self.rawComponent("G", z);
                F = self.rawComponent("F", z);
            elseif self.evp.primaryComponent == "F"
                F = self.rawComponent("F", z);
            end

            switch normalization
                case Normalization.uMax
                    factors = max(abs(F), [], 1);
                case Normalization.wMax
                    factors = max(abs(G), [], 1);
                case Normalization.geostrophic
                    if isempty(G)
                        error("InternalModesBasisSet:UnsupportedNormalization", ...
                            "Geostrophic normalization requires a G component.");
                    end
                    weight = self.solver.N2(z)/self.solver.g;
                    factors = sqrt(trapz(z, weight(:).*G.*G, 1));
                case Normalization.omegaConstant
                    if isempty(F)
                        error("InternalModesBasisSet:UnsupportedNormalization", ...
                            "omegaConstant normalization requires an F component.");
                    end
                    factors = sqrt(trapz(z, F.*F, 1)/diff(self.solver.zDomain));
                case Normalization.kConstant
                    if isempty(G)
                        factors = sqrt(trapz(z, F.*F, 1));
                    else
                        weight = max(0,self.solver.N2(z) - self.solver.f0*self.solver.f0)/self.solver.g;
                        factors = sqrt(trapz(z, weight(:).*G.*G, 1));
                    end
                otherwise
                    error("InternalModesBasisSet:UnsupportedNormalization", ...
                        "Unsupported normalization requested.");
            end
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
            gram = self.partialGramMatrix(component, self.solver.zDomain(1), self.solver.zDomain(2));
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
            self.requireSolver("partialGramMatrix");
            gram = self.solver.componentGramMatrix(self, component, [zMin zMax]);
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
                self InternalModesBasisSet
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
                self InternalModesBasisSet
                coefficientsA (:,1) double
                coefficientsB (:,1) double
                options.component {mustBeTextScalar} = self.evp.primaryComponent
            end

            gram = self.gramMatrix(options.component);
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
        end
    end

    methods (Access = private)
        function values = rawComponent(self, component, z)
            self.requireSolver("evaluate");
            component = string(component);
            switch component
                case self.evp.primaryComponent
                    values = self.solver.evaluateNativeModes(self.nativeModes, z);
                case "F"
                    if self.evp.primaryComponent ~= "G"
                        self.unsupported("evaluate F diagnostic");
                    end
                    derivative = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
                    values = derivative .* self.h;
                otherwise
                    self.unsupported("evaluate " + component);
            end
        end

        function requireSolver(self, operationName)
            if isempty(self.solver)
                error("InternalModesBasisSet:UnsupportedOperation", ...
                    "The basis set does not have a solver reference and cannot perform %s.", operationName);
            end
        end

        function unsupported(~, operationName)
            error("InternalModesBasisSet:UnsupportedOperation", ...
                "InternalModesBasisSet does not support %s for this basis.", operationName);
        end
    end

    methods (Static)
        function basisSet = constantStratification(options)
            % Create an analytical constant-stratification basis placeholder.
            %
            % Unsupported general basis-set operations throw structured
            % errors until exact analytical implementations are added.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = InternalModesBasisSet.constantStratification(options)
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: analytical basis placeholder
            arguments
                options.metadata struct = struct()
            end

            metadata = options.metadata;
            metadata.analyticalBasis = "constantStratification";
            basisSet = InternalModesBasisSet(metadata=metadata);
        end

        function basisSet = exponentialStratification(options)
            % Create an analytical exponential-stratification basis placeholder.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = InternalModesBasisSet.exponentialStratification(options)
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: analytical basis placeholder
            arguments
                options.metadata struct = struct()
            end

            metadata = options.metadata;
            metadata.analyticalBasis = "exponentialStratification";
            basisSet = InternalModesBasisSet(metadata=metadata);
        end

        function basisSet = wkbApproximation(options)
            % Create an analytical WKB basis placeholder.
            %
            % - Topic: Create basis sets
            % - Declaration: basisSet = InternalModesBasisSet.wkbApproximation(options)
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: analytical basis placeholder
            arguments
                options.metadata struct = struct()
            end

            metadata = options.metadata;
            metadata.analyticalBasis = "wkbApproximation";
            basisSet = InternalModesBasisSet(metadata=metadata);
        end
    end
end
