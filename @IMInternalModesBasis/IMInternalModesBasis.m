classdef IMInternalModesBasis < IMBasisSet
    % Store solved internal-mode basis functions.
    %
    % `IMInternalModesBasis` evaluates the physical variables with explicit
    % `F(z)` and `G(z)` methods. If the EVP solves `G`, then `F(z)` is
    % recovered by `evp.FfromGz`, whose default relation is
    % $$F_j(z)=h_j\frac{\partial G_j}{\partial z}(z).$$
    % If the EVP solves `F`, then `G(z)` is recovered by `evp.GfromFz`.
    % The default inverse is hydrostatic, but wave factories install
    % wave-specific inverse relations.
    % Normalization is shared across both variables: if a rule gives scale
    % $$s_j$$, then both diagnostic variables for mode $$j$$ are divided by
    % that same factor. Standard rules are installed on the basis set, and
    % custom rules are added after solving with `addNormalization`.
    %
    % ```matlab
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet.normalization = Normalization.geostrophic;
    % F = basisSet.F(z);
    % G = basisSet.G(z);
    % ```
    %
    % - Topic: Create internal-mode bases
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

        % Buoyancy frequency squared profile.
        %
        % `N2` has signature `values = N2(z)`. It is copied from the
        % internal-mode EVP unless supplied explicitly by an analytical
        % basis class.
        %
        % - Topic: Inspect internal-mode bases
        N2
    end

    methods
        function self = IMInternalModesBasis(options)
            % Create an internal-mode basis set.
            %
            % - Topic: Create internal-mode bases
            % - Declaration: basisSet = IMInternalModesBasis(options)
            % - Parameter options.solver: solver reference
            % - Parameter options.evp: internal-mode EVP descriptor
            % - Parameter options.nativeModes: native mode columns
            % - Parameter options.eigenvalues: retained eigenvalues
            % - Parameter options.h: equivalent depths
            % - Parameter options.modeNumber: physical mode numbers
            % - Parameter options.modeSelectionDiagnostics: mode-selection diagnostics
            % - Parameter options.normalization: active normalization rule name or enum value
            % - Parameter options.metadata: additional metadata
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.N2: buoyancy frequency squared function
            % - Returns basisSet: internal-mode basis set
            arguments
                options.solver = []
                options.evp IMInternalModes
                options.nativeModes double = zeros(0,0)
                options.eigenvalues double = zeros(1,0)
                options.h double = zeros(1,0)
                options.modeNumber double = []
                options.modeSelectionDiagnostics struct = struct()
                options.normalization = []
                options.metadata struct = struct()
                options.zDomain (1,2) double = [NaN NaN]
                options.N2 = []
            end

            self@IMBasisSet(solver=options.solver, evp=options.evp, nativeModes=options.nativeModes, eigenvalues=options.eigenvalues, modeNumber=options.modeNumber, modeSelectionDiagnostics=options.modeSelectionDiagnostics, normalization=options.normalization, metadata=options.metadata, zDomain=options.zDomain);
            if isempty(options.N2)
                self.N2 = options.evp.N2;
            else
                self.N2 = options.N2;
            end
            if ~isa(self.N2, "function_handle")
                error("IMInternalModesBasis:InvalidN2", "N2 must be a function handle.");
            end

            self.h = reshape(options.h,1,[]);
            if isempty(self.h) && ~isempty(self.eigenvalues)
                self.h = options.evp.hFromEigenvalue(reshape(self.eigenvalues,1,[]));
                self.h = reshape(self.h,1,[]);
            end
            if ~isempty(self.h)
                nModes = max([size(self.nativeModes,2), length(self.eigenvalues), length(self.modeNumber)]);
                if length(self.h) ~= nModes
                    error("IMInternalModesBasis:InvalidEquivalentDepthCount", "h must contain one equivalent depth for each retained mode.");
                end
            end

            self = self.installInternalModeNormalizationRules();
            if isempty(options.normalization)
                if self.evp.modeFamily == "geostrophic"
                    self.normalization = "geostrophic";
                elseif string(self.evp.name) == "waveModesAtWavenumber"
                    self.normalization = "kConstant";
                else
                    self.normalization = "unity";
                end
            end
        end

        function G = G(self, z, options)
            % Evaluate `G` modes.
            %
            % If the EVP formulation is `G`, this evaluates the solved
            % canonical variable. If the formulation is `F`, `G` is recovered
            % by `evp.GfromFz`. The default relation is
            % $$G_j(z)=-\frac{g}{N^2(z)}
            % \frac{\partial F_j}{\partial z}(z),$$
            % but individual EVPs may supply a different diagnostic
            % relation.
            %
            % - Topic: Evaluate internal-mode bases
            % - Declaration: G = G(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule name or enum value
            % - Returns G: evaluated `G` modes
            arguments
                self IMInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            G = self.rawVariable("G", z) ./ self.normalizationFactors(options.normalization);
        end

        function F = F(self, z, options)
            % Evaluate `F` modes.
            %
            % If the EVP formulation is `F`, this evaluates the solved
            % canonical variable. If the formulation is `G`, `F` is recovered
            % by `evp.FfromGz`; the default relation is
            % $$F_j(z)=h_j\frac{\partial G_j}{\partial z}(z).$$
            %
            % - Topic: Evaluate internal-mode bases
            % - Declaration: F = F(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule name or enum value
            % - Returns F: evaluated `F` modes
            arguments
                self IMInternalModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end

            F = self.rawVariable("F", z) ./ self.normalizationFactors(options.normalization);
        end

        function gram = gramMatrix(self, options)
            % Return a Gram matrix for `F` or `G`.
            %
            % With no arguments this uses the solved formulation over the full
            % basis-set domain. Use `variable="F"` or `variable="G"` to choose
            % a physical variable, and `zBounds=[zMin zMax]` to restrict the
            % interior integral. Endpoint terms are included only when the
            % interval contains the corresponding physical endpoint:
            % $$M_{ij}=\int_{z_a}^{z_b} w(z)V_i(z)V_j(z)\,dz+
            % \text{included endpoint terms}.$$
            % The requested variable must have a known inner product; if it
            % does not, this method throws
            % `IMInternalModesBasis:UnavailableInnerProduct` rather than
            % returning an incomplete Gram matrix.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: gram = gramMatrix(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns gram: Gram matrix
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            self.validateZBounds(options.zBounds(1), options.zBounds(2));
            gram = self.variableGramMatrix(string(options.variable), options.zBounds, true);
        end

        function windowModes = partialWindowModes(self, options)
            % Diagonalize a partial-depth Gram matrix for `F` or `G`.
            %
            % This computes the eigendecomposition of the symmetric Gram
            % matrix on `zBounds`. If `variable` is omitted, the solved
            % formulation is used. The requested variable must have a known
            % inner product.
            %
            % - Topic: Analyze Gram matrices
            % - Declaration: windowModes = partialWindowModes(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns windowModes: window-mode decomposition
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            self.validateZBounds(options.zBounds(1), options.zBounds(2));
            gram = self.variableGramMatrix(string(options.variable), options.zBounds, true);
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
            % used. The requested variable must have a known inner product.
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
            % used. The requested variable must have a known inner product.
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

            gram = self.gramMatrix(variable=options.variable);
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
            FReferences = self.modeSignReferences("F", zSurface, zGrid);
            GReferences = self.modeSignReferences("G", zSurface, zGrid);
            signs = ones(1,size(self.nativeModes,2));
            for iMode = 1:size(self.nativeModes,2)
                reference = FReferences(iMode);
                if ~isfinite(reference) || reference == 0
                    reference = GReferences(iMode);
                end
                if isfinite(reference) && reference < 0
                    signs(iMode) = -1;
                end
            end
            self.nativeModes = self.nativeModes .* signs;
        end

        function factor = innerProductNormFactor(self, iMode, options)
            % Return the `F` or `G` inner-product norm factor.
            %
            % This is the raw factor
            % $$s_j=\sqrt{|\langle V_j,V_j\rangle|}$$ for `variable` equal
            % to `F` or `G`. If `variable` is omitted, the solved
            % formulation is used. The requested variable must have a known
            % inner product. Custom normalization rules
            % registered with `addNormalization` call this method.
            %
            % - Topic: Developer topics
            % - Declaration: factor = innerProductNormFactor(basisSet,iMode,options)
            % - Parameter iMode: retained mode index
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns factor: raw inner-product norm factor
            % - Developer: true
            arguments
                self IMInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            variable = string(options.variable);
            gram = self.variableGramMatrix(variable, self.zDomain, false);
            factor = sqrt(abs(gram(iMode,iMode)));
        end

        function factor = geostrophicNormFactor(self, iMode)
            % Return the hydrostatic geostrophic normalization factor.
            %
            % This rule is installed for `modeFamily="geostrophic"` and
            % chooses one shared scale factor for each coupled `F`/`G`
            % mode. Baroclinic modes use
            % $$s_j^2=\langle G_j,G_j\rangle_G,$$
            % so normalized modes satisfy
            % $$\langle G_j,G_j\rangle_G=1,\qquad
            % \langle F_j,F_j\rangle_F=h_j.$$
            % A barotropic zero mode uses the `F` norm divided by
            % $$\sqrt{z_\mathrm{surface}-z_\mathrm{bottom}}$$ and is a
            % separate null-mode convention.
            %
            % - Topic: Developer topics
            % - Declaration: factor = geostrophicNormFactor(basisSet,iMode)
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
            % This is $$s_j=\max_z |V_j^{\mathrm{raw}}(z)|$$ for the
            % requested variable on the basis-set integration grid. If
            % `variable` is omitted, the solved formulation is used.
            %
            % - Topic: Developer topics
            % - Declaration: factor = maxAbsFactor(basisSet,iMode,options)
            % - Parameter iMode: retained mode index
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns factor: maximum absolute variable amplitude
            % - Developer: true
            arguments
                self IMInternalModesBasis
                iMode (1,1) double {mustBeInteger, mustBePositive}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            variable = string(options.variable);
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

    methods (Static)
        function basisSet = constantStratification(options)
            % Create an analytical constant-stratification basis set.
            %
            % - Topic: Create internal-mode bases
            % - Declaration: basisSet = IMInternalModesBasis.constantStratification(options)
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

            basisSet = IMBasisSetConstantStratification(evp=options.evp, N0=options.N0, zDomain=options.zDomain, nModes=options.nModes, normalization=options.normalization, metadata=options.metadata);
        end

        function basisSet = exponentialStratification(options)
            % Create an analytical exponential-stratification basis set.
            %
            % - Topic: Create internal-mode bases
            % - Declaration: basisSet = IMInternalModesBasis.exponentialStratification(options)
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

            basisSet = IMBasisSetExponentialStratification(evp=options.evp, N0=options.N0, b=options.b, zDomain=options.zDomain, nModes=options.nModes, normalization=options.normalization, metadata=options.metadata);
        end
    end

    methods (Access = protected)
        function context = context(self)
            context = context@IMBasisSet(self);
            context.N2 = @(z) self.N2(z);
        end

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
            self.unsupported("evaluate analytical derivatives without a subclass implementation");
        end

        function values = rawVariable(self, variable, z)
            arguments
                self IMInternalModesBasis
                variable {mustBeTextScalar, mustBeMember(variable, ["F", "G"])}
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            variable = string(variable);
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
                    dGdz = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
                    values = self.evp.FfromGz(z(:), dGdz, self.h, self.context());
                case "F"
                    if variable ~= "G"
                        self.unsupported("evaluate " + variable);
                    end
                    dFdz = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
                    values = self.evp.GfromFz(z(:), dFdz, self.h, self.context());
                otherwise
                    self.unsupported("evaluate " + variable);
            end
        end

        function gram = variableGramMatrix(self, variable, zBounds, useNormalized)
            z = self.innerProductGrid(zBounds);
            spec = self.evp.innerProduct(variable);
            IMInternalModesBasis.assertInnerProductAvailable(spec);
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
                    integrand = weight(:).*values(:,iMode).*values(:,jMode);
                    value = self.integrateInnerProduct(z, integrand, zBounds);
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
    end

    methods (Access = private)
        function self = installInternalModeNormalizationRules(self)
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

        function references = modeSignReferences(self, variable, zSurface, zGrid)
            surfaceValues = [];
            gridValues = [];
            try
                surfaceValues = self.rawVariable(variable, zSurface);
                gridValues = self.rawVariable(variable, zGrid);
            catch exception
                if string(exception.identifier) ~= "IMBasisSet:UnsupportedOperation"
                    rethrow(exception)
                end
            end

            references = NaN(1,size(self.nativeModes,2));
            for iMode = 1:size(self.nativeModes,2)
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
                references(iMode) = value;
            end
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
    end

    methods (Static, Access = private)
        function assertInnerProductAvailable(spec)
            if isfield(spec, "hasInnerProduct") && spec.hasInnerProduct
                return;
            end
            error("IMInternalModesBasis:UnavailableInnerProduct", "The %s inner product is unavailable for this EVP and cannot be used as a Gram matrix. %s", string(spec.variable), string(spec.reason));
        end

    end

end
