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
    % - Topic: Evaluate modes
    % - Topic: Analyze modes
    % - Topic: Inspect basis sets
    % - Topic: Developer topics
    % - Declaration: classdef IMInternalModesBasis < IMBasisSet

    properties (SetAccess = protected)
        % Equivalent depths for the retained internal modes.
        %
        % These are computed from the parent EVP as
        % $$h_j=\texttt{evp.hFromEigenvalue}(\lambda_j),$$
        % where $$\lambda_j$$ is `eigenvalues(j)`.
        %
        % - Topic: Inspect basis sets
        h

        % Buoyancy frequency squared profile.
        %
        % `N2` has signature `values = N2(z)` and is copied from the
        % internal-mode EVP.
        %
        % - Topic: Inspect basis sets
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
            % - Parameter options.modeNumber: physical mode numbers
            % - Parameter options.modeSelectionDiagnostics: mode-selection diagnostics
            % - Parameter options.normalization: active normalization rule name or enum value
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: internal-mode basis set
            arguments
                options.solver IMSolver
                options.evp IMInternalModes
                options.nativeModes (:,:) double
                options.eigenvalues (1,:) double {mustBeReal, mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.modeSelectionDiagnostics struct = struct()
                options.normalization = []
                options.metadata struct = struct()
            end

            self@IMBasisSet(solver=options.solver, evp=options.evp, nativeModes=options.nativeModes, eigenvalues=options.eigenvalues, modeNumber=options.modeNumber, modeSelectionDiagnostics=options.modeSelectionDiagnostics, normalization=options.normalization, metadata=options.metadata);
            self.N2 = options.evp.N2;
            self.h = reshape(options.evp.hFromEigenvalue(self.eigenvalues),1,[]);
            if length(self.h) ~= size(self.nativeModes,2)
                error("IMInternalModesBasis:InvalidEquivalentDepthCount", "hFromEigenvalue must return one equivalent depth for each native mode column.");
            end

            self = self.addNormalization("uMax", @(basisSet,iMode) basisSet.maxAmplitudeNormFactor(iMode, variable="F"));
            self = self.addNormalization("wMax", @(basisSet,iMode) basisSet.maxAmplitudeNormFactor(iMode, variable="G"));
            self = self.addNormalization("surfacePressure", @(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode));
            if self.evp.modeFamily == "hydrostatic"
                self = self.addNormalization("geostrophic", @(basisSet,iMode) basisSet.geostrophicNormFactor(iMode));
            end
            if string(self.evp.name) == "waveModesAtWavenumber"
                self = self.addNormalization("kConstant", @(basisSet,iMode) basisSet.innerProductNormFactor(iMode, variable="G"));
            end
            if isempty(options.normalization)
                if self.evp.modeFamily == "hydrostatic"
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
            % - Topic: Evaluate modes
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
            % - Topic: Evaluate modes
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
            %
            % $$
            % M_{ij}=\int_{z_a}^{z_b} w(z)V_i(z)V_j(z)\,dz+
            % \sum_\ell \gamma_\ell L_\ell[V_i]L_\ell[V_j]+
            % \sum_\ell \alpha_\ell V_i(z_\ell)V_j(z_\ell).
            % $$
            %
            % Use `endpointGramTerms` to inspect the prepared endpoint
            % vectors that generate the rank-one endpoint updates.
            % The requested variable must have a known inner product; if it
            % does not, this method throws
            % `IMInternalModesBasis:UnavailableInnerProduct` rather than
            % returning an incomplete Gram matrix.
            %
            % - Topic: Analyze modes
            % - Declaration: gram = gramMatrix(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns gram: Gram matrix
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end
            gram = self.variableGramMatrix(string(options.variable), options.zBounds, true);
        end

        function terms = endpointGramTerms(self, options)
            % Prepare rank-one endpoint terms for `F` or `G` Gram matrices.
            %
            % `endpointGramTerms` returns endpoint vectors over all retained
            % modes. Solved-form endpoint weights use
            %
            % $$
            % L_\ell[V_j]=c_\ell V_j(z_\ell)-d_\ell p(z_\ell)\frac{\partial V_j}{\partial z}(z_\ell),
            % $$
            %
            % and contribute
            %
            % $$
            % M \leftarrow M+\gamma_\ell L_\ell L_\ell^\mathsf{T}.
            % $$
            %
            % Catalog endpoint value terms use $$V_j(z_\ell)$$ and
            % contribute
            %
            % $$
            % M \leftarrow M+\alpha_\ell V_\ell V_\ell^\mathsf{T}.
            % $$
            %
            % Endpoint terms are omitted when `zBounds` does not include
            % that endpoint.
            %
            % - Topic: Developer topics — Gram-matrix assembly
            % - Declaration: terms = endpointGramTerms(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Parameter options.normalization: normalization rule name or enum value
            % - Parameter options.useNormalized: whether returned values use the normalization
            % - Returns terms: struct array with `location`, `coefficient`, `values`, and `kind`
            % - Developer: true
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
                options.normalization = self.normalization
                options.useNormalized (1,1) logical = true
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end

            variable = string(options.variable);
            spec = self.evp.innerProduct(variable);
            if ~isfield(spec, "hasInnerProduct") || ~spec.hasInnerProduct
                error("IMInternalModesBasis:UnavailableInnerProduct", "The %s inner product is unavailable for this EVP and cannot be used as a Gram matrix. %s", string(spec.variable), string(spec.reason));
            end
            terms = struct("location", {}, "coefficient", {}, "values", {}, "kind", {});
            if variable == self.evp.formulation
                terms = endpointGramTerms@IMBasisSet(self, zBounds=options.zBounds, normalization=options.normalization, useNormalized=options.useNormalized);
            end

            normalizationFactors = [];
            for iTerm = 1:numel(spec.endpointInnerProductTerms)
                endpointTerm = spec.endpointInnerProductTerms(iTerm);
                if string(endpointTerm.location) == "surface"
                    zEndpoint = self.zDomain(2);
                else
                    zEndpoint = self.zDomain(1);
                end

                tolerance = 100*eps(max(1,max(abs([options.zBounds(:); zEndpoint]))));
                endpointIsIncluded = abs(min(options.zBounds) - zEndpoint) <= tolerance || abs(max(options.zBounds) - zEndpoint) <= tolerance;
                if ~endpointIsIncluded
                    continue;
                end

                values = self.rawVariable(endpointTerm.variable, zEndpoint);
                if options.useNormalized
                    if isempty(normalizationFactors)
                        normalizationFactors = self.normalizationFactors(options.normalization);
                    end
                    values = values ./ normalizationFactors;
                end
                terms(end+1) = struct("location", string(endpointTerm.location), "coefficient", endpointTerm.coefficient, "values", values, "kind", "endpointInnerProductTerm");
            end
        end

        function windowModes = partialWindowModes(self, options)
            % Diagonalize a partial-depth Gram matrix for `F` or `G`.
            %
            % This computes the eigendecomposition of the symmetric Gram
            % matrix on `zBounds`. If `variable` is omitted, the solved
            % formulation is used. The requested variable must have a known
            % inner product.
            %
            % - Topic: Analyze modes
            % - Declaration: windowModes = partialWindowModes(basisSet,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Parameter options.zBounds: integration bounds `[zMin zMax]`
            % - Returns windowModes: window-mode decomposition
            arguments
                self IMInternalModesBasis
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
                options.zBounds (1,2) double {mustBeReal, mustBeFinite} = self.zDomain
            end

            if options.zBounds(1) >= options.zBounds(2)
                error("IMBasisSet:InvalidInterval", "zBounds must be increasing.");
            end
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
            % - Topic: Analyze modes
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
            % - Topic: Analyze modes
            % - Declaration: spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB,options)
            % - Parameter coefficientsA: first modal coefficients
            % - Parameter coefficientsB: second modal coefficients
            % - Parameter options.variable: optional variable name, `"F"` or `"G"`
            % - Returns spectrum: modal cross-spectrum
            arguments
                self IMInternalModesBasis
                coefficientsA (:,1) double {mustBeFinite}
                coefficientsB (:,1) double {mustBeFinite}
                options.variable {mustBeTextScalar, mustBeMember(options.variable, ["F", "G"])} = self.evp.formulation
            end

            nModes = size(self.nativeModes,2);
            if length(coefficientsA) ~= nModes || length(coefficientsB) ~= nModes
                error("IMBasisSet:InvalidCoefficientCount", "Coefficient vectors must contain one value for each retained mode (%d).", nModes);
            end
            gram = self.gramMatrix(variable=options.variable);
            spectrum = diag(gram).*real(coefficientsA(:).*conj(coefficientsB(:)));
        end
    end

    methods
        function self = orientModeSigns(self)
            % Orient modes so the surface `F` value is positive when possible.
            %
            % This developer utility applies the internal-mode sign
            % convention used after numerical solves: prefer a finite
            % nonzero surface `F` value, fall back to the largest `F`
            % value on the solver grid, then fall back to `G`. The same
            % sign flip is applied to the coupled `F`/`G` mode pair.
            % `IMInternalModesBasis` is a value class, so callers must keep
            % the returned basis set:
            %
            % ```matlab
            % basisSet = basisSet.orientModeSigns();
            % ```
            %
            % - Topic: Developer topics — Diagnostic variables
            % - Declaration: basisSet = orientModeSigns(basisSet)
            % - Returns basisSet: basis set with oriented native mode signs
            % - Developer: true
            if isempty(self.nativeModes)
                return;
            end
            zSurface = self.zDomain(2);
            zGrid = self.solver.innerProductGrid(self.zDomain);
            variables = ["F", "G"];
            references = NaN(2,size(self.nativeModes,2));
            for iVariable = 1:numel(variables)
                surfaceValues = [];
                gridValues = [];
                try
                    surfaceValues = self.rawVariable(variables(iVariable), zSurface);
                    gridValues = self.rawVariable(variables(iVariable), zGrid);
                catch exception
                    if string(exception.identifier) ~= "IMBasisSet:UnsupportedOperation"
                        rethrow(exception)
                    end
                end

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
                    references(iVariable,iMode) = value;
                end
            end
            signs = ones(1,size(self.nativeModes,2));
            for iMode = 1:size(self.nativeModes,2)
                reference = references(1,iMode);
                if ~isfinite(reference) || reference == 0
                    reference = references(2,iMode);
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
            % The factor is computed from raw, unnormalized modes before
            % the active basis normalization is applied.
            %
            % - Topic: Developer topics — Normalization rules
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
            % This developer utility implements the normalization rule
            % installed for `modeFamily="hydrostatic"`. It chooses one
            % shared raw scale factor for each coupled `F`/`G` mode.
            % Baroclinic modes use
            % $$s_j^2=\langle G_j,G_j\rangle_G,$$
            % so normalized modes satisfy
            % $$\langle G_j,G_j\rangle_G=1,\qquad
            % \langle F_j,F_j\rangle_F=h_j.$$
            % A barotropic zero mode uses the `F` norm divided by
            % $$\sqrt{z_\mathrm{surface}-z_\mathrm{bottom}}$$ and is a
            % separate null-mode convention.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = geostrophicNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: raw geostrophic scale factor
            % - Developer: true
            if self.modeNumber(iMode) == 0
                gram = self.variableGramMatrix("F", self.zDomain, false);
                factor = sqrt(abs(gram(iMode,iMode)))/sqrt(diff(self.zDomain));
                return;
            end
            factor = self.innerProductNormFactor(iMode, variable="G");
        end

        function factor = maxAmplitudeNormFactor(self, iMode, options)
            % Return the maximum amplitude of `F` or `G`.
            %
            % This is $$s_j=\max_z |V_j^{\mathrm{raw}}(z)|$$ for the
            % requested variable on the basis-set integration grid. If
            % `variable` is omitted, the solved formulation is used.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = maxAmplitudeNormFactor(basisSet,iMode,options)
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
            z = self.solver.innerProductGrid(self.zDomain);
            values = self.rawVariable(variable, z);
            factor = max(abs(values(:,iMode)));
        end

        function factor = surfacePressureNormFactor(self, iMode)
            % Return the raw surface `F` value.
            %
            % This developer utility returns the raw surface value
            % $$s_j=F_j^{\mathrm{raw}}(z_\mathrm{surface}).$$
            % It gives unit surface `F` value when the raw surface value is
            % finite and nonzero. If that value is unavailable, zero, or
            % nonfinite, it returns `1` so normalization remains well
            % defined.
            %
            % - Topic: Developer topics — Normalization rules
            % - Declaration: factor = surfacePressureNormFactor(basisSet,iMode)
            % - Parameter iMode: retained mode index
            % - Returns factor: raw surface-pressure scale factor
            % - Developer: true
            values = self.rawVariable("F", self.zDomain(2));
            factor = values(1,iMode);
            if ~isfinite(factor) || abs(factor) <= 1e-12
                factor = 1;
            end
        end

        function values = rawVariable(self, variable, z)
            % Evaluate raw physical `F` or `G` modes.
            %
            % This developer utility evaluates unnormalized internal-mode
            % physical variables. The public `F` and `G` methods apply the
            % active basis normalization after this step:
            % $$F_j(z)=F_j^{\mathrm{raw}}(z)/s_j,\qquad
            % G_j(z)=G_j^{\mathrm{raw}}(z)/s_j.$$
            % If `variable` is the solved formulation, values come
            % directly from the native modes. Otherwise the diagnostic
            % variable is recovered through the EVP-owned diagnostic
            % relation.
            %
            % - Topic: Developer topics — Diagnostic variables
            % - Declaration: values = rawVariable(basisSet,variable,z)
            % - Parameter variable: `"F"` or `"G"`
            % - Parameter z: physical coordinate
            % - Returns values: unnormalized physical mode values
            % - Developer: true
            arguments
                self IMInternalModesBasis
                variable {mustBeTextScalar, mustBeMember(variable, ["F", "G"])}
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            variable = string(variable);
            context = self.evp.contextForSolver(self.solver);
            if variable == self.evp.formulation
                values = self.solver.evaluateNativeModes(self.nativeModes, z);
                return;
            end

            switch self.evp.formulation
                case "G"
                    if variable ~= "F"
                        self.unsupported("evaluate " + variable);
                    end
                    dGdz = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
                    values = self.evp.FfromGz(z(:), dGdz, self.h, context);
                case "F"
                    if variable ~= "G"
                        self.unsupported("evaluate " + variable);
                    end
                    dFdz = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
                    values = self.evp.GfromFz(z(:), dFdz, self.h, context);
                otherwise
                    self.unsupported("evaluate " + variable);
            end
        end
    end

    methods (Access = protected)
        function gram = variableGramMatrix(self, variable, zBounds, useNormalized)
            z = self.solver.innerProductGrid(zBounds);
            context = self.evp.contextForSolver(self.solver);
            spec = self.evp.innerProduct(variable);
            if ~isfield(spec, "hasInnerProduct") || ~spec.hasInnerProduct
                error("IMInternalModesBasis:UnavailableInnerProduct", "The %s inner product is unavailable for this EVP and cannot be used as a Gram matrix. %s", string(spec.variable), string(spec.reason));
            end
            normalizationFactors = [];
            if useNormalized
                normalizationFactors = self.normalizationFactors(self.normalization);
                values = self.rawVariable(variable, z) ./ normalizationFactors;
            else
                values = self.rawVariable(variable, z);
            end
            weight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, context);
            if isscalar(weight)
                weight = weight*ones(size(z));
            end
            gram = zeros(size(values,2), size(values,2));
            for iMode = 1:size(values,2)
                for jMode = iMode:size(values,2)
                    integrand = weight(:).*values(:,iMode).*values(:,jMode);
                    value = self.solver.integrateInnerProduct(z, integrand, zBounds);
                    gram(iMode,jMode) = value;
                    gram(jMode,iMode) = value;
                end
            end

            endpointTerms = self.endpointGramTerms(variable=variable, zBounds=zBounds, useNormalized=useNormalized);
            for iTerm = 1:numel(endpointTerms)
                valuesAtEndpoint = endpointTerms(iTerm).values(:);
                gram = gram + endpointTerms(iTerm).coefficient*(valuesAtEndpoint*valuesAtEndpoint.');
            end
        end
    end

end
