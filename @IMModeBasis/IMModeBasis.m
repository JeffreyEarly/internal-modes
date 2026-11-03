classdef IMModeBasis
    % Evaluate one mode family with shared identities and normalization names.
    %
    % Factories in `IMModes` create these immutable continuous bases. Wave
    % profiles obey $$F_j = D G_j'$$, where $$D = z_s-z_b$$. Every
    % normalization multiplies all fields in one column by the same
    % positive factor; signed metric signatures remain explicit.
    % Scientific and construction records preserve the defining problem,
    % selected columns and backend evidence. Shared capabilities describe
    % available operations and their requirements; observation-grid rank
    % and recovery evidence belong to a prepared transform. Supported
    % scalar physical recipes can bind this family to caller-supplied
    % samples through projectionOnGrid or scalarTransform; continuous
    % construction remains independent of sampling and recovery policy.
    % `IMVerticalSampling` evaluates prescribed point/cell observations
    % without changing the basis or choosing numerical recovery.
    % Zero-APV bases also expose physical coefficient forms and immutable
    % rotations within the retained continuous subspace. Canonical endpoint
    % labels identify source rows; the source-coordinate map identifies
    % their current normalized or rotated linear combinations.
    %
    % ```matlab
    % basis = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=8);
    % F = basis.F(z);
    % ```
    %
    % - Topic: Inspect modes
    % - Topic: Evaluate modes
    % - Topic: Project modes
    % - Topic: Form zero-APV quadratic forms
    % - Topic: Rotate zero-APV modes
    % - Topic: Developer topics
    % - Declaration: classdef IMModeBasis

    properties (SetAccess = private)
        % Scientific family name.
        % - Topic: Inspect modes
        family (1,1) string
        % Complete retained column count, including special modes.
        % - Topic: Inspect modes
        nColumns (1,1) double
        % Stable scientific column identities.
        % - Topic: Inspect modes
        columnLabels (1,:) string
        % One row per column, with branch and branch-local index.
        % - Topic: Inspect modes
        modeInfo table
        % Resolved normalization name.
        % - Topic: Inspect modes
        normalization (1,1) string
        % Authoritative evaluation, normalization and operation capabilities.
        % - Topic: Inspect modes
        capabilities (1,1) struct
        % Scientific family, pairing, boundaries and selected column identity.
        % - Topic: Inspect modes
        scientificSpecification (1,1) struct
        % Requested and resolved backend settings and construction evidence.
        % - Topic: Inspect modes
        constructionRecord (1,1) struct
        % Immutable physical column configuration.
        % - Topic: Inspect modes
        background IMModeBackground
        % Increasing bottom and surface coordinates, in meters.
        % - Topic: Inspect modes
        zDomain (1,2) double
        % Signs remaining after signed normalization, otherwise positive one.
        % - Topic: Inspect modes
        normSignatures (1,:) double
        % Horizontal identity: fixedWavenumber, perColumn, independent or unspecified.
        % - Topic: Inspect modes
        horizontalKind (1,1) string
        % Fixed nonnegative horizontal wavenumber, empty for other identities.
        % - Topic: Inspect modes
        k (1,:) double
        % Signed horizontal roots in inverse square meters, otherwise empty.
        % - Topic: Inspect modes
        kSquared (1,:) double
        % Physical modal frequencies when defined by the family.
        % - Topic: Inspect modes
        frequency (1,:) double
        % Caller extension metadata, separate from canonical scientific records.
        % - Topic: Inspect modes
        metadata (1,1) struct
        % Canonical endpoint row identities for the zero-APV source.
        % - Topic: Inspect modes
        canonicalEndpointLabels (1,:) string = strings(1,0)
        % Rotation creation evidence, separate from ordinary EVP spectra.
        %
        % The parent-coordinate pencil and map describe rotation creation.
        % Selection narrows aligned arrays; positive profile normalization
        % leaves this snapshot unchanged. Current forms and sourceCoordinateMap
        % describe the resulting coefficient coordinates.
        % - Topic: Rotate zero-APV modes
        rotationMetadata (1,1) struct
    end

    properties (Dependent, SetAccess = private)
        % Equivalent depths; unavailable for families without equivalent depth.
        % - Topic: Inspect modes
        h
        % EVP eigenvalues; unavailable for families without an EVP.
        % - Topic: Inspect modes
        eigenvalues
        % Zero-APV-only map from canonical endpoints to normalized columns.
        %
        % The nCanonicalEndpoints-by-nColumns map $$M$$ obeys
        % $$A_{\mathrm{current}}=A_{\mathrm{canonical}}M$$ for every profile,
        % derivative and endpoint trace. Canonical row identities remain
        % separate from mixed column labels. Other families throw when this
        % unavailable coordinate map is requested.
        % - Topic: Inspect modes
        sourceCoordinateMap
        % Physical-energy form in current zero-APV coefficient coordinates.
        %
        % For coefficients $$a$$, the restricted energy self-pairing is
        % $$a^* H a$$. Under a coordinate change $$T$$, the form becomes
        % $$H_{\mathrm{new}}=T^* H T$$. The state convention includes
        % the existing zero-APV factors of k, f0 and g.
        % - Topic: Form zero-APV quadratic forms
        energyMatrix
        % Surface-buoyancy form in current zero-APV coordinates.
        % - Topic: Form zero-APV quadratic forms
        surfaceBuoyancyMatrix
        % Bottom-buoyancy form in current zero-APV coordinates.
        % - Topic: Form zero-APV quadratic forms
        bottomBuoyancyMatrix
        % Canonical endpoint-response form in current coordinates.
        % - Topic: Form zero-APV quadratic forms
        endpointResponseMetric
    end

    properties (Access = private)
        source
        columns (1,:) double
        includesBarotropic (1,1) logical = false
        waveConvention (1,1) logical
        intrinsicEvaluator
        equivalentDepths_ (1,:) double
        eigenvalues_ (1,:) double
        hasEquivalentDepths_ (1,1) logical = false
        hasEigenvalues_ (1,1) logical = false
        hasFrequency_ (1,1) logical = false
        profileVariables_ (1,:) string
        % Immutable derived scales depend on source directions and normalization.
        scaleFactors_ (1,:) double
        sourceUnityFactors_ (1,:) double
        integrationIntervals_ (:,2) double
        % Original source inventory fixes the normalization integration budget.
        normalizationColumnCount_ (1,1) double
        zeroAPVColumnMap_ (:,:) double = zeros(0)
        zeroAPVSourceRotation_ (:,:) double = zeros(0)
        zeroAPVRightGram_ (:,:) double = zeros(0)
        zeroAPVEnergyBoundary_ (1,1) struct = struct()
    end

    methods
        function self = IMModeBasis(options)
            % Wrap selected continuous columns in the common family contract.
            %
            % Numerical and analytical internal-mode sources are evaluated
            % with their `unity` normalization. The resulting column factors
            % are computed once. Immutable construction inputs make these
            % cached factors valid for every subsequent evaluation.
            %
            % - Topic: Developer topics
            % - Parameter options.source: existing numerical or analytical basis with the same declared stratification as background, or empty
            % - Parameter options.background: common physical column configuration
            % - Parameter options.family: scientific family name
            % - Parameter options.columns: selected source columns; all by default
            % - Parameter options.modeInfo: table with columnLabel, branch and branchIndex
            % - Parameter options.normalization: resolved common normalization name
            % - Parameter options.includeBarotropic: prepend the analytical rigid-lid F0=1, G0=0 sector to a G source
            % - Parameter options.waveConvention: whether to enforce F = D G prime
            % - Parameter options.horizontalKind: declared horizontal identity; unspecified for an unassigned direct wrapper
            % - Parameter options.k: nonnegative scalar for fixed-wavenumber identity, otherwise empty
            % - Parameter options.kSquared: signed per-column horizontal roots for fixed-frequency identity
            % - Parameter options.h: optional equivalent depths for selected columns
            % - Parameter options.eigenvalues: optional selected eigenvalues
            % - Parameter options.hasEquivalentDepths: advertise an empty equivalent-depth spectrum
            % - Parameter options.hasEigenvalues: advertise an empty EVP spectrum
            % - Parameter options.frequency: optional physical modal frequencies
            % - Parameter options.hasFrequency: advertise an empty physical-frequency spectrum
            % - Parameter options.intrinsicEvaluator: function of variable, z and derivative order
            % - Parameter options.capabilities: intrinsic variable/derivative declarations such as F, Fz, Gz and Gzz
            % - Parameter options.variables: optional intrinsic variable declarations; G also declares Gz/Gzz
            % - Parameter options.metadata: opaque caller extension metadata
            % - Parameter options.integrationIntervals: exact intrinsic support intervals
            % - Parameter options.scientificInput: factory-owned scientific request and resolved EVP declaration
            % - Parameter options.constructionInput: factory-owned backend request, resolved settings and selection evidence
            % - Returns self: immutable continuous family basis
            % - Developer: true
            arguments (Input)
                options.source = []
                options.background (1,1) IMModeBackground
                options.family (1,1) string
                options.columns (1,:) double {mustBeInteger,mustBePositive} = []
                options.modeInfo table = table()
                options.normalization (1,1) string = "depthRMSF"
                options.includeBarotropic (1,1) logical = false
                options.waveConvention (1,1) logical = false
                options.horizontalKind (1,1) string {mustBeMember(options.horizontalKind,["fixedWavenumber","perColumn","independent","unspecified"])} = "unspecified"
                options.k (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
                options.kSquared (1,:) double {mustBeReal,mustBeFinite} = []
                options.h (1,:) double {mustBeReal} = []
                options.eigenvalues (1,:) double {mustBeReal} = []
                options.hasEquivalentDepths (1,1) logical = false
                options.hasEigenvalues (1,1) logical = false
                options.frequency (1,:) double {mustBeReal} = []
                options.hasFrequency (1,1) logical = false
                options.intrinsicEvaluator = []
                options.capabilities (1,:) string = ["F","G","Gz","Gzz"]
                options.variables (1,:) string = strings(1,0)
                options.metadata (1,1) struct = struct()
                options.integrationIntervals (:,2) double {mustBeReal,mustBeFinite} = zeros(0,2)
                options.scientificInput (1,1) struct = struct()
                options.constructionInput (1,1) struct = struct()
            end
            arguments (Output)
                self (1,1) IMModeBasis
            end
            self.source = options.source;
            self.background = options.background;
            self.family = options.family;
            self.zDomain = options.background.zDomain;
            self.normalization = options.normalization;
            self.includesBarotropic = options.includeBarotropic;
            self.waveConvention = options.waveConvention;
            self.k = options.k;
            self.kSquared = options.kSquared;
            self.horizontalKind = options.horizontalKind;
            if self.horizontalKind == "unspecified"
                if ~isempty(self.k)
                    self.horizontalKind = "fixedWavenumber";
                elseif ~isempty(self.kSquared)
                    self.horizontalKind = "perColumn";
                end
            end
            self.frequency = options.frequency;
            self.hasFrequency_ = options.hasFrequency || ~isempty(options.frequency);
            self.metadata = options.metadata;
            self.rotationMetadata = struct(name="",isRotated=false,eigenvalues=zeros(1,0),rightMetricSignatures=zeros(1,0),leftMetricSignatures=zeros(1,0),normalizationConvention="",residuals=struct(),h0=zeros(1,0),parentColumnLabels=strings(1,0),coordinateChange=zeros(0),leftMatrix=zeros(0),rightMatrix=zeros(0));
            self.integrationIntervals_ = options.integrationIntervals;
            if ~isempty(self.integrationIntervals_) && (any(self.integrationIntervals_(:,1) >= self.integrationIntervals_(:,2)) || any(self.integrationIntervals_(:,1) < self.zDomain(1)) || any(self.integrationIntervals_(:,2) > self.zDomain(2)) || any(self.integrationIntervals_(2:end,1) < self.integrationIntervals_(1:end-1,2)))
                error("IMModeBasis:InvalidIntegrationIntervals","Integration intervals must be increasing, disjoint and contained in zDomain.");
            end
            self.intrinsicEvaluator = options.intrinsicEvaluator;
            if ~isempty(self.intrinsicEvaluator) && ~isa(self.intrinsicEvaluator,"function_handle")
                error("IMModeBasis:InvalidEvaluator","intrinsicEvaluator must be a function handle.");
            end
            if isempty(self.source)
                nSource = height(options.modeInfo);
                self.profileVariables_ = options.capabilities;
                if ~isempty(options.variables)
                    self.profileVariables_ = options.variables;
                    if ismember("G",self.profileVariables_)
                        self.profileVariables_ = unique([self.profileVariables_,"Gz","Gzz"],"stable");
                    end
                end
                if nSource > 0 && isempty(self.intrinsicEvaluator)
                    error("IMModeBasis:MissingEvaluator","Nonempty intrinsic modes require an intrinsicEvaluator.");
                end
            elseif isa(self.source,"IMNeutralLabelModesBasis")
                nSource = numel(self.source.columnLabels);
                self.profileVariables_ = ["G","Gz","Gzz","APV"];
            elseif isa(self.source,"IMInternalModesBasis") || isa(self.source,"IMAnalyticalInternalModesBasis")
                nSource = numel(self.source.eigenvalues);
                self.profileVariables_ = ["F","G","Gz","Gzz"];
                self.sourceUnityFactors_ = self.source.normalizationFactors("unity");
                factors = self.sourceUnityFactors_;
                self.source = self.source.addNormalization("unity",@(~,j) factors(j));
                self.source.normalization = "unity";
            elseif isa(self.source,"IMGeostrophicZeroAPVModesBasis") || isa(self.source,"IMAnalyticalGeostrophicZeroAPVModesBasis")
                if numel(self.source.k) ~= 1
                    error("IMModeBasis:MultipleWavenumbers","Wrap one zero-APV wavenumber page at a time.");
                end
                nSource = numel(self.source.endpoints);
                self.profileVariables_ = ["F","G","Gz","Gzz"];
            else
                error("IMModeBasis:UnsupportedSource","source must be a supported continuous mode basis.");
            end
            if ~isempty(self.source) && ~isequal(self.source.zDomain,self.zDomain)
                error("IMModeBasis:IncompatibleDomain","The source and background must have the same zDomain.");
            end
            if self.includesBarotropic
                if isempty(self.source) || self.source.evp.formulation ~= "G" || self.family ~= "hydrostaticF" || self.background.surfaceCondition ~= "rigidLid"
                    error("IMModeBasis:InvalidBarotropicExtension","An analytical barotropic extension requires a rigid-lid hydrostatic G source.");
                end
                nSource = nSource+1;
            end
            self.normalizationColumnCount_ = nSource;
            self.columns = options.columns;
            if isempty(self.columns)
                self.columns = 1:nSource;
            elseif any(self.columns > nSource) || numel(unique(self.columns)) ~= numel(self.columns)
                error("IMModeBasis:InvalidColumns","Select distinct existing source columns.");
            end
            self.nColumns = numel(self.columns);
            if self.family == "geostrophicZeroAPV"
                self = self.initializeZeroAPVCoordinates(options.scientificInput);
            end
            self.hasEquivalentDepths_ = options.hasEquivalentDepths;
            self.hasEigenvalues_ = options.hasEigenvalues;
            self.modeInfo = options.modeInfo;
            required = ["columnLabel","branch","branchIndex"];
            if width(self.modeInfo) == 0 && self.nColumns == 0
                self.modeInfo = table(strings(0,1),strings(0,1),zeros(0,1),VariableNames=required);
            end
            if height(self.modeInfo) ~= self.nColumns || ~all(ismember(required,string(self.modeInfo.Properties.VariableNames)))
                error("IMModeBasis:InvalidModeInfo","modeInfo must have one row per selected column and columnLabel, branch and branchIndex variables.");
            end
            self.columnLabels = reshape(string(self.modeInfo.columnLabel),1,[]);
            if numel(unique(self.columnLabels)) ~= self.nColumns
                error("IMModeBasis:InvalidColumnLabels","Column labels must be distinct scientific identities.");
            end
            if ~isempty(options.h)
                self.equivalentDepths_ = options.h;
                self.hasEquivalentDepths_ = true;
            elseif ~isempty(self.source) && isprop(self.source,"h")
                spectrum = self.source.h;
                if self.includesBarotropic, spectrum = [Inf spectrum]; end
                self.equivalentDepths_ = spectrum(self.columns);
                self.hasEquivalentDepths_ = true;
            end
            if ~isempty(options.eigenvalues)
                self.eigenvalues_ = options.eigenvalues;
                self.hasEigenvalues_ = true;
            elseif ~isempty(self.source) && isprop(self.source,"eigenvalues")
                spectrum = self.source.eigenvalues;
                if self.includesBarotropic, spectrum = [0 spectrum]; end
                self.eigenvalues_ = spectrum(self.columns);
                self.hasEigenvalues_ = true;
            end
            if (self.hasEquivalentDepths_ && numel(self.equivalentDepths_) ~= self.nColumns) || (self.hasEigenvalues_ && numel(self.eigenvalues_) ~= self.nColumns) || (~isempty(self.frequency) && numel(self.frequency) ~= self.nColumns)
                error("IMModeBasis:InvalidSpectrumCount","h, eigenvalues and frequency must have one value per selected column when supplied.");
            end
            if self.hasFrequency_ && numel(self.frequency) ~= self.nColumns
                error("IMModeBasis:InvalidSpectrumCount","A declared physical-frequency spectrum requires one value per retained column.");
            end
            if self.waveConvention && ~isempty(self.source) && self.source.evp.formulation ~= "G"
                error("IMModeBasis:WaveFormulation","The standard wave convention requires a G-form source.");
            end
            if self.horizontalKind == "fixedWavenumber"
                if ~isscalar(self.k) || ~isempty(self.kSquared)
                    error("IMModeBasis:InvalidHorizontalIdentity","A fixed-wavenumber basis requires one k and no per-column kSquared.");
                end
            elseif self.horizontalKind == "perColumn"
                if ~isempty(self.k) || numel(self.kSquared) ~= self.nColumns
                    error("IMModeBasis:InvalidHorizontalIdentity","Per-column dispersion requires one signed kSquared per retained column and no scalar k.");
                end
            elseif ~isempty(self.k) || ~isempty(self.kSquared)
                error("IMModeBasis:InvalidHorizontalIdentity","Independent or unspecified horizontal identity cannot carry a solved k or kSquared.");
            end
            self.capabilities = self.buildCapabilities(options.scientificInput);
            [self.scaleFactors_,self.normSignatures] = self.normalizationScales();
            [self.scientificSpecification,self.constructionRecord] = self.createRecords(options.scientificInput,options.constructionInput);
            if self.family == "geostrophicZeroAPV", self = self.refreshZeroAPVRecords(); end
        end

        function values = F(self,z)
            % Evaluate F profiles as a sample-by-column matrix.
            % - Topic: Evaluate modes
            % - Parameter z: physical coordinate column
            % - Returns values: normalized F profiles
            values = self.evaluate(z,"F");
        end

        function values = G(self,z)
            % Evaluate G profiles as a sample-by-column matrix.
            % - Topic: Evaluate modes
            % - Parameter z: physical coordinate column
            % - Returns values: normalized G profiles
            values = self.evaluate(z,"G");
        end

        function values = Gz(self,z)
            % Evaluate physical first derivatives of G.
            % - Topic: Evaluate modes
            % - Parameter z: physical coordinate column
            % - Returns values: normalized first derivatives
            values = self.evaluate(z,"G",1);
        end

        function values = Gzz(self,z)
            % Evaluate classical physical second derivatives of G.
            % - Topic: Evaluate modes
            % - Parameter z: physical coordinate column
            % - Returns values: normalized second derivatives
            values = self.evaluate(z,"G",2);
        end

        function values = APV(self,z)
            % Evaluate the family-defined APV channel when available.
            % - Topic: Evaluate modes
            % - Parameter z: physical coordinate column
            % - Returns values: normalized APV profiles
            values = self.evaluate(z,"APV");
        end

        function values = evaluate(self,z,variable,derivativeOrder)
            % Evaluate a profile or its physical vertical derivatives.
            %
            % Values have one row per point and one column per retained
            % mode. Neutral-interval derivatives use inside endpoint traces.
            %
            % - Topic: Evaluate modes
            % - Parameter z: physical coordinate column within zDomain
            % - Parameter variable: F, Fz, G, Gz, Gzz or APV, when available
            % - Parameter derivativeOrder: physical derivative order, zero by default
            % - Returns values: normalized sample-by-column matrix
            arguments (Input)
                self (1,1) IMModeBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                variable (1,1) string {mustBeMember(variable,["F","Fz","G","Gz","Gzz","APV"])}
                derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
            end
            arguments (Output)
                values (:,:) double
            end
            for alias = self.capabilities.aliases
                if variable == alias.variable
                    variable = alias.baseVariable;
                    derivativeOrder = derivativeOrder+alias.derivativeOffset;
                    break
                end
            end
            variableIndex = find(self.capabilities.variables == variable,1);
            if isempty(variableIndex)
                error("IMModeBasis:UnavailableVariable","The %s family has no %s profile.",self.family,variable);
            end
            if ~ismember(derivativeOrder,self.capabilities.derivativeOrders{variableIndex})
                error("IMModeBasis:UnsupportedDerivative","The %s family does not support derivative order %d of %s.",self.family,derivativeOrder,variable);
            end
            if any(z < self.zDomain(1)) || any(z > self.zDomain(2))
                error("IMModeBasis:OutsideDomain","Evaluation points must lie within zDomain.");
            end
            values = self.unscaledValues(variable,z,derivativeOrder)./self.scaleFactors_;
        end

        function value = get.h(self)
            if ~self.hasEquivalentDepths_
                error("IMModeBasis:UnavailableEquivalentDepth","The %s family has no equivalent depth.",self.family);
            end
            value = self.equivalentDepths_;
        end

        function value = get.eigenvalues(self)
            if ~self.hasEigenvalues_
                error("IMModeBasis:UnavailableEigenvalues","The %s family is not defined by an EVP.",self.family);
            end
            value = self.eigenvalues_;
        end

        function map = get.sourceCoordinateMap(self)
            if self.family ~= "geostrophicZeroAPV"
                error("IMModeBasis:UnavailableCoordinateMap","A canonical endpoint coordinate map is defined only for the geostrophicZeroAPV family.");
            end
            map = self.zeroAPVColumnMap_./self.scaleFactors_;
            if any(~isfinite(map),"all")
                error("IMModeBasis:UnrepresentableCoordinates","The normalized canonical coordinate map is not finite. Choose a representable normalization or pencil.");
            end
        end

        function matrix = get.energyMatrix(self)
            matrix = self.zeroAPVForm("energyMatrix");
        end

        function matrix = get.surfaceBuoyancyMatrix(self)
            matrix = self.zeroAPVForm("surfaceBuoyancyMatrix");
        end

        function matrix = get.bottomBuoyancyMatrix(self)
            matrix = self.zeroAPVForm("bottomBuoyancyMatrix");
        end

        function matrix = get.endpointResponseMetric(self)
            matrix = self.zeroAPVForm("endpointResponseMetric");
        end

        function projection = projectionOnGrid(self,z,weights,options)
            % Construct a family projection on fixed physical samples.
            %
            % Supported scalar F/G channels retain their full continuous
            % signed forms and literal source-observation requirements.
            % Neutral-label APV uses its prescribed continuous dual for the
            % fixed-endpoint sine ladder. Physical projection applies its
            % continuous-normalized dual; sample recovery is independent.
            %
            % - Topic: Project modes
            % - Parameter z: increasing physical sample column
            % - Parameter weights: fixed physical dz quadrature weights
            % - Parameter options.variable: supported F/G channel, or APV for neutral labels
            % - Parameter options.columns: selected facade columns, all by default
            % - Parameter options.prepareRecovery: prepare the legacy physical-metric recovery, or return only the physical record
            % - Parameter options.failureMode: diagnostic or error for unresolved recovery
            % - Parameter options.rankTolerance: dimensionless relative recovery rank cutoff
            % - Returns projection: family map from sampled values to normalized-basis coefficients
            arguments (Input)
                self (1,1) IMModeBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                weights (:,1) double {mustBeReal,mustBeFinite}
                options.variable (1,1) string = ""
                options.columns (1,:) double {mustBeInteger,mustBePositive} = 1:self.nColumns
                options.prepareRecovery (1,1) logical = true
                options.failureMode (1,1) string {mustBeMember(options.failureMode,["diagnostic","error"])} = "diagnostic"
                options.rankTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-12
            end
            arguments (Output)
                projection (1,1) IMProjection
            end
            selectedColumns = options.columns;
            if any(selectedColumns > self.nColumns) || numel(unique(selectedColumns)) ~= numel(selectedColumns)
                error("IMModeBasis:InvalidColumns","Select distinct existing basis columns.");
            end
            if ~isa(self.source,"IMNeutralLabelModesBasis")
                recipe = self.projectionRecipe(variable=options.variable);
                projection = IMProjection.fromRecipe(recipe,z,weights,columns=selectedColumns,prepareRecovery=options.prepareRecovery,failureMode=options.failureMode,rankTolerance=options.rankTolerance);
                return
            end
            if ~ismember(options.variable,["","APV"])
                error("IMModeBasis:UnavailableProjection","Neutral labels provide the prescribed physical APV dual only.");
            end
            if isempty(selectedColumns)
                % Validate the physical rule through an empty source, without
                % preparing discarded APV columns or their physical norms.
                source = IMNeutralLabelModesBasis(self.source.neutralIntervals,zDomain=self.zDomain,nModes=0,rotationRate=self.source.rotationRate,latitude=self.source.latitude,g0=self.source.g0,gd=self.source.gd);
            else
                source = self.source;
            end
            base = source.projectionOnGrid(z,weights,columns=self.columns(selectedColumns));
            factors = self.scaleFactors_(selectedColumns);
            forward = factors(:).*base.forwardMatrix;
            majorant = base.majorantGramMatrix./(factors(:)*factors);
            provenance = base.provenance;
            provenance.normalization = self.normalization;
            provenance.continuousAPVGramMatrix = majorant;
            projection = IMProjection.fromPrescribedDual(forward,eye(numel(selectedColumns)),eye(numel(selectedColumns)),majorantGramMatrix=majorant,columnLabels=self.columnLabels(selectedColumns),provenance=provenance);
        end
    end

    methods
        basis = selectColumns(self,columns)
        basis = withNormalization(self,name)
        matrix = generalizedEnergyMatrix(self,options)
        basis = rotateBoundaryDepth(self,options)
        basis = rotateSurfaceBuoyancy(self,options)
        basis = rotateWithPencil(self,options)
        transform = discreteTransform(self,grid)
        transform = scalarTransform(self,z,weights,options)
        [z,design] = modeRootPoints(self,generatingMode)
    end

    methods (Hidden)
        [HF,HG] = hydrostaticNorms(self)
        info = samplingIntegrationInfo(self)
        function recipe = projectionRecipe(self,options)
            % Describe the family projection metric and continuous reference.
            % - Topic: Project modes
            % - Parameter options.variable: supported scalar F/G channel, or APV for neutral labels
            % - Returns recipe: family-owned reference and capability metadata
            % - Developer: true
            arguments (Input)
                self (1,1) IMModeBasis
                options.variable (1,1) string = ""
            end
            arguments (Output)
                recipe (1,1) struct
            end
            variable = options.variable;
            if isa(self.source,"IMNeutralLabelModesBasis") && ismember(variable,["","APV"])
                recipe = self.source.projectionRecipe(variable="APV");
                factors = self.scaleFactors_;
                recipe.columnLabels = self.columnLabels;
                recipe.normalization = self.normalization;
                recipe.targetGramMatrix = eye(self.nColumns);
                recipe.majorantGramMatrix = recipe.majorantGramMatrix(self.columns,self.columns)./(factors(:)*factors);
                recipe.activeColumnMask = true(1,self.nColumns);
                recipe.evaluateFunction = @(z) self.APV(z);
                recipe.provenance.normalization = self.normalization;
                recipe.provenance.continuousAPVGramMatrix = recipe.majorantGramMatrix;
                return
            end
            supported = ismember(self.family,["hydrostaticF","hydrostaticG","geostrophicAPV","meanDensityAnomaly"]);
            if ~supported
                error("IMModeBasis:UnavailableProjection","The %s family has no common scalar physical projection in this adapter. Evaluable profiles may still support an explicitly prepared observation recovery.",self.family);
            end
            if variable == ""
                if self.family == "hydrostaticF", variable = "F"; else, variable = "G"; end
            end
            if ~ismember(variable,["F","G"]) || (self.family == "meanDensityAnomaly" && variable ~= "G")
                error("IMModeBasis:UnavailableProjection","The %s family does not define an independent scalar physical %s pairing.",self.family,variable);
            end
            active = true(1,self.nColumns);
            hasSource = ~isempty(self.source);
            if hasSource
                self.assertValueOnlyPairing(self.source.evp.innerProduct(variable),variable);
                recipe = self.source.projectionRecipe(variable=variable);
                if ~recipe.isAvailable
                    error("IMModeBasis:UnavailableProjection","%s",recipe.reason);
                end
                target = recipe.targetGramMatrix;
                majorant = recipe.majorantGramMatrix;
                active = recipe.activeColumnMask;
                if self.includesBarotropic
                    if variable == "F", constantNorm = self.background.D; else, constantNorm = 0; end
                    target = blkdiag(constantNorm,target);
                    majorant = blkdiag(constantNorm,majorant);
                    active = [variable == "F",active];
                    recipe.provenance.barotropicAugmentation = "F0=1, G0=0; exact Dirichlet G traces imply zero F0/internal cross pairing";
                end
                target = target(self.columns,self.columns);
                majorant = majorant(self.columns,self.columns);
                active = active(self.columns);
            else
                pairings = self.scientificSpecification.resolved.pairings;
                if ~isfield(pairings,variable) || ~pairings.(variable).hasInnerProduct
                    error("IMModeBasis:UnavailableProjection","The %s family has no retained declaration of the requested %s pairing.",self.family,variable);
                end
                spec = pairings.(variable);
                if ~isfield(spec,"surfaceWeights"), spec.surfaceWeights = struct([]); end
                if ~isfield(spec,"bottomWeights"), spec.bottomWeights = struct([]); end
                if ~isfield(spec,"endpointInnerProductTerms"), spec.endpointInnerProductTerms = struct([]); end
                if variable == "F", weightFunction = @(z) ones(size(z)); else, weightFunction = @(z) self.background.N2(z)/self.background.g; end
                target = zeros(self.nColumns);
                if self.family == "hydrostaticF" && self.background.surfaceCondition == "rigidLid"
                    if variable == "F", target = diag(self.background.D*double(self.modeInfo.branch == "barotropic")); end
                elseif self.nColumns > 0
                    error("IMModeBasis:UnavailableProjection","A nonempty scalar family requires an existing continuous source form.");
                end
                majorant = target;
                provenance = struct(representation="intrinsic",referenceIntegration="empty form or exact analytical barotropic norm",solveAccuracy="unverified",referenceConvergence="unverified");
                recipe = IMProjection.createRecipe(variable=variable,columnLabels=self.columnLabels,targetGramMatrix=target,majorantGramMatrix=majorant,provenance=provenance,zDomain=self.zDomain,weightFunction=weightFunction,spec=spec);
            end
            self.assertValueOnlyPairing(recipe.spec,variable);
            if self.family == "hydrostaticF" && self.background.surfaceCondition == "rigidLid" && variable == "G"
                active = self.modeInfo.branch.' ~= "barotropic";
            end
            factors = self.scaleFactors_;
            recipe.columnLabels = self.columnLabels;
            recipe.normalization = self.normalization;
            recipe.targetGramMatrix = target./(factors(:)*factors);
            recipe.majorantGramMatrix = majorant./(factors(:)*factors);
            recipe.activeColumnMask = active;
            recipe.evaluateFunction = @(z) self.evaluate(z,variable);
            recipe.provenance.family = self.family;
            recipe.provenance.projection = "sampledGalerkin";
            recipe.provenance.normalization = self.normalization;
            recipe.provenance.columnLabels = self.columnLabels;
            recipe.provenance.sourceColumns = self.columns;
            recipe.provenance.targetCoordinateChange = "selected source unity profiles divided by common normalization scales";
        end
    end

    methods (Access = private)
        capabilities = buildCapabilities(self,scientificInput)
        [scientific,construction] = createRecords(self,scientificInput,constructionInput)
        self = initializeZeroAPVCoordinates(self,scientificInput)
        matrix = zeroAPVForm(self,name,isRaw)
        boundary = zeroAPVEnergyBoundary(self,g0,gd)
        self = refreshZeroAPVRecords(self)
        basis = rotateZeroAPVPencil(self,name,left,right,ordering,convention,g0,gd)
        bound = zeroAPVNormBound(self,boundary)

        function assertValueOnlyPairing(~,spec,variable)
            if ~spec.hasInnerProduct
                error("IMModeBasis:UnavailableProjection","%s",spec.reason);
            end
            terms = [spec.surfaceWeights;spec.bottomWeights];
            if ~isempty(terms) && any([terms.coefficient] ~= 0 & [terms.d] ~= 0)
                error("IMModeBasis:UnavailableProjection","The %s pairing requires independent derivative endpoint observations; scalar values cannot represent it.",variable);
            end
            terms = spec.endpointInnerProductTerms;
            if ~isempty(terms) && any([terms.coefficient] ~= 0 & string({terms.variable}) ~= variable)
                error("IMModeBasis:UnavailableProjection","The %s pairing requires independent companion-variable endpoint observations; scalar values cannot represent it.",variable);
            end
        end

        function values = unscaledValues(self,variable,z,derivativeOrder)
            if self.nColumns == 0
                values = zeros(numel(z),0);
                return
            end
            if self.family == "geostrophicZeroAPV"
                if variable == "F" && derivativeOrder == 0
                    raw = self.source.F(z);
                elseif variable == "F" && derivativeOrder == 1
                    raw = -self.source.N2(z)/self.background.g.*self.source.G(z);
                elseif variable == "G" && derivativeOrder == 0
                    raw = self.source.G(z);
                elseif variable == "G" && derivativeOrder == 1
                    raw = -self.background.g*self.k^2/self.background.f0^2*self.source.F(z);
                elseif variable == "G" && derivativeOrder == 2
                    raw = self.k^2*self.source.N2(z)/self.background.f0^2.*self.source.G(z);
                else
                    error("IMModeBasis:UnsupportedDerivative","Zero-APV profiles support F through first and G through second derivatives.");
                end
                values = (raw/self.zeroAPVSourceRotation_)*self.zeroAPVColumnMap_;
                return
            end
            if isempty(self.source)
                values = self.intrinsicEvaluator(variable,z,derivativeOrder);
            elseif isa(self.source,"IMNeutralLabelModesBasis")
                if variable == "G" && derivativeOrder <= 2
                    names = ["G","Gz","Gzz"];
                    values = self.source.(names(derivativeOrder+1))(z);
                elseif variable == "APV" && derivativeOrder == 0
                    values = self.source.APV(z);
                else
                    error("IMModeBasis:UnsupportedDerivative","This neutral-label derivative is unavailable.");
                end
            elseif self.waveConvention && variable == "F"
                if derivativeOrder > 1
                    error("IMModeBasis:UnsupportedDerivative","Wave F supports derivative orders zero and one.");
                end
                values = diff(self.zDomain)*self.sourceGDerivative(z,derivativeOrder+1);
            elseif variable == "G" && derivativeOrder > 0
                values = self.sourceGDerivative(z,derivativeOrder);
            elseif derivativeOrder == 0
                if isempty(self.sourceUnityFactors_)
                    values = self.source.(variable)(z);
                else
                    values = self.source.(variable)(z,normalization="unity");
                end
            elseif variable == "F" && derivativeOrder == 1
                if isprop(self.source,"evp") && self.source.evp.formulation == "F"
                    values = self.source.uz(z,normalization="unity");
                else
                    values = -self.background.N2(z)/self.background.g.*self.unscaledSourceG(z);
                end
            else
                error("IMModeBasis:UnsupportedDerivative","This profile derivative is unavailable.");
            end
            if self.includesBarotropic
                constant = zeros(numel(z),1);
                if variable == "F" && derivativeOrder == 0, constant(:) = 1; end
                values = [constant values];
            end
            values = values(:,self.columns);
            if self.includesBarotropic && variable == "G" && derivativeOrder == 0
                values(z == self.zDomain(1) | z == self.zDomain(2),:) = 0;
            end
        end

        function values = unscaledSourceG(self,z)
            if isempty(self.sourceUnityFactors_)
                values = self.source.G(z);
            else
                values = self.source.G(z,normalization="unity");
            end
        end

        function values = sourceGDerivative(self,z,order)
            if order > 2
                error("IMModeBasis:UnsupportedDerivative","G supports derivative orders zero through two.");
            end
            if order == 1 && ismethod(self.source,"Gz")
                if isempty(self.sourceUnityFactors_)
                    values = self.source.Gz(z);
                else
                    values = self.source.Gz(z,normalization="unity");
                end
            elseif isa(self.source,"IMBasisSet") && self.source.evp.formulation == "G"
                values = self.source.solver.evaluatePhysicalDerivative(self.source.nativeModes,z,order)./self.sourceUnityFactors_;
            elseif isprop(self.source,"evp") && self.source.evp.formulation == "G" && order == 1
                values = self.source.uz(z,normalization="unity");
            elseif isprop(self.source,"evp") && self.source.evp.formulation == "G"
                evp = self.source.evp;
                if evp.name == "waveModesAtWavenumber"
                    potential = evp.parameters.k^2;
                    weight = (self.background.N2(z)-self.background.f0^2)/self.background.g;
                elseif evp.name == "waveModesAtFrequency"
                    potential = 0;
                    weight = (self.background.N2(z)-evp.parameters.omega^2)/self.background.g;
                else
                    potential = 0;
                    weight = self.background.N2(z)/self.background.g;
                end
                values = (potential-weight.*self.source.eigenvalues).*self.unscaledSourceG(z);
            elseif isprop(self.source,"evp")
                if order == 1
                    values = self.source.F(z,normalization="unity")./self.source.h;
                else
                    values = -self.background.N2(z)/self.background.g.*self.unscaledSourceG(z)./self.source.h;
                end
            else
                % Zero-APV endpoint responses obey the regular first-order pair.
                factor = -self.background.g*self.k^2/self.background.f0^2;
                if order == 1
                    values = factor*self.source.F(z);
                else
                    values = -factor*self.background.N2(z)/self.background.g.*self.source.G(z);
                end
            end
        end

        function [scales,signatures] = normalizationScales(self)
            names = ["depthRMSF","depthRMSG","waveEnergy","generalizedEnergy","generalizedPotentialEnstrophy","endpointResponse","maxAbsF","maxAbsG","rotationRightMetric"];
            if ~ismember(self.normalization,names)
                error("IMModeBasis:UnknownNormalization","Choose one of the common normalization names: %s.",join(names,", "));
            end
            if ~ismember(self.normalization,self.capabilities.normalizationNames)
                error("IMModeBasis:UnavailableNormalization","The %s family does not provide %s with its declared profiles and horizontal identity.",self.family,self.normalization);
            end
            signatures = ones(1,self.nColumns);
            if self.nColumns == 0
                scales = zeros(1,0);
                return
            end
            switch self.normalization
                case "endpointResponse"
                    if ~(isa(self.source,"IMGeostrophicZeroAPVModesBasis") || isa(self.source,"IMAnalyticalGeostrophicZeroAPVModesBasis"))
                        error("IMModeBasis:UnavailableNormalization","endpointResponse is defined only for canonical zero-APV endpoint responses.");
                    end
                    scales = ones(1,self.nColumns);
                    return
                case "rotationRightMetric"
                    squares = reshape(diag(self.zeroAPVRightGram_),1,[]);
                    scales = sqrt(abs(squares));
                    signatures = sign(squares);
                case "generalizedEnergy"
                    if self.family == "geostrophicZeroAPV"
                        boundary = self.zeroAPVEnergyBoundary([],[]);
                        gram = self.zeroAPVForm("energyMatrix",true);
                        majorant = gram;
                        if isfinite(boundary.g0), gram = gram+boundary.g0*self.zeroAPVForm("surfaceBuoyancyMatrix",true); end
                        if isfinite(boundary.gd), gram = gram+boundary.gd*self.zeroAPVForm("bottomBuoyancyMatrix",true); end
                        if isfinite(boundary.g0), majorant = majorant+abs(boundary.g0)*self.zeroAPVForm("surfaceBuoyancyMatrix",true); end
                        if isfinite(boundary.gd), majorant = majorant+abs(boundary.gd)*self.zeroAPVForm("bottomBuoyancyMatrix",true); end
                        squares = reshape(diag(gram),1,[]);
                        majorantSquares = reshape(diag(majorant),1,[]);
                        arithmeticBound = self.zeroAPVNormBound(boundary);
                        relativeNullTolerance = 1e-12;
                        if any(~isfinite(majorantSquares)) || any(~isfinite(arithmeticBound)) || any(majorantSquares <= 0) || any(abs(squares) <= relativeNullTolerance*max(majorantSquares,arithmeticBound))
                            error("IMModeBasis:NullNormalization","A retained zero-APV column has an unresolved signed-energy norm relative to its physical majorant or canonical arithmetic cancellation bound (relative tolerance %g). Select a non-null span or another normalization.",relativeNullTolerance);
                        end
                        signatures = sign(squares);
                        scales = sqrt(abs(squares));
                    else
                        if isempty(self.sourceUnityFactors_) || ~ismethod(self.source,"gramMatrix") || self.waveConvention
                            error("IMModeBasis:UnavailableNormalization","This family has no generalized-energy G pairing.");
                        end
                        gram = self.source.gramMatrix(variable="G");
                        squares = reshape(diag(gram),1,[]);
                        if self.includesBarotropic, squares = [0 squares]; end
                        squares = squares(self.columns);
                        signatures = sign(squares);
                        scales = sqrt(abs(squares));
                    end
                case "generalizedPotentialEnstrophy"
                    if isempty(self.sourceUnityFactors_) || ~ismember("generalizedPotentialEnstrophy",self.source.normalizationNames())
                        error("IMModeBasis:UnavailableNormalization","This family has no generalized-potential-enstrophy pairing.");
                    end
                    factors = self.source.normalizationFactors("generalizedPotentialEnstrophy")./self.sourceUnityFactors_;
                    scales = factors(self.columns);
                otherwise
                    if isa(self.source,"IMNeutralLabelModesBasis") && self.normalization == "depthRMSG"
                        scales = ones(1,self.nColumns);
                        return
                    end
                    if isa(self.source,"IMNeutralLabelModesBasis") && self.normalization == "maxAbsG"
                        lengths = diff(self.source.neutralIntervals,1,2);
                        scales = sqrt(2*diff(self.zDomain)./reshape(lengths(self.source.intervalIndex(self.columns)),1,[]));
                        return
                    end
                    [z,weights] = self.normalizationQuadrature();
                    D = diff(self.zDomain);
                    switch self.normalization
                        case {"depthRMSF","maxAbsF"}
                            self.requireVariable("F");
                            values = self.unscaledValues("F",z,0);
                        case {"depthRMSG","maxAbsG"}
                            self.requireVariable("G");
                            values = self.unscaledValues("G",z,0);
                        case "waveEnergy"
                            if ~self.waveConvention
                                error("IMModeBasis:UnavailableNormalization","waveEnergy requires the standard wave F/G convention.");
                            end
                            F = self.unscaledValues("F",z,0);
                            G = self.unscaledValues("G",z,0);
                            scales = sqrt(sum(weights.*(abs(F).^2+(self.k*D)^2*abs(G).^2),1)/D);
                    end
                    if startsWith(self.normalization,"depthRMS")
                        scales = sqrt(sum(weights.*abs(values).^2,1)/D);
                    elseif startsWith(self.normalization,"maxAbs")
                        scales = self.maximumAmplitudes(z,values);
                    end
            end
            if any(~isfinite(scales)) || any(scales <= 0)
                error("IMModeBasis:NullNormalization","A retained mode has zero or nonfinite %s magnitude. Choose an admissible normalization; modes are never silently removed.",self.normalization);
            end
        end

        function scales = maximumAmplitudes(self,z,values)
            variable = extractAfter(self.normalization,"maxAbs");
            scales = max(abs(values),[],1);
            options = optimset("TolX",max(1e-12,eps(diff(self.zDomain))),"Display","off");
            for j = 1:self.nColumns
                magnitude = abs(values(:,j));
                candidates = find(magnitude(2:end-1) >= magnitude(1:end-2) & magnitude(2:end-1) >= magnitude(3:end) & (magnitude(2:end-1) > magnitude(1:end-2) | magnitude(2:end-1) > magnitude(3:end)))+1;
                for i = reshape(candidates,1,[])
                    [~,minimum] = fminbnd(@(zz) -self.columnMagnitude(variable,zz,j),z(i-1),z(i+1),options);
                    scales(j) = max(scales(j),-minimum);
                end
            end
        end

        function value = columnMagnitude(self,variable,z,column)
            values = self.unscaledValues(variable,z,0);
            value = abs(values(1,column));
        end

        function requireVariable(self,variable)
            if ~ismember(variable,self.capabilities.variables)
                error("IMModeBasis:UnavailableNormalization","%s requires a %s profile, which this family does not provide.",self.normalization,variable);
            end
        end

        function [z,weights] = normalizationQuadrature(self,includeExtrema)
            if nargin < 2, includeExtrema = true; end
            if ~isempty(self.source) && isprop(self.source,"solver") && isa(self.source.solver,"IMSolverSpectral")
                solver = self.source.solver;
                bounds = [min(solver.xNative),max(solver.xNative)];
                [x,nativeWeights] = legpts(max(128,2*solver.nEVP),bounds);
                z = solver.zOfX(x);
                % The degree-one Chebyshev representation of x gives its
                % exact physical derivative dx/dz for the solver's map.
                coordinateCoefficients = zeros(solver.nEVP,1);
                coordinateCoefficients(1:2) = [mean(bounds);diff(bounds)/2];
                jacobian = solver.evaluatePhysicalDerivative(coordinateCoefficients,z,1);
                weights = nativeWeights(:)./jacobian;
                if includeExtrema && startsWith(self.normalization,"maxAbs")
                    z = unique([z;self.zDomain.']);
                end
                return
            end
            if ~isempty(self.source) && isprop(self.source,"solver")
                [z,weights] = self.source.solver.nativeQuadratureRule(self.zDomain);
                return
            end
            if ~isempty(self.integrationIntervals_)
                bounds = self.integrationIntervals_;
            elseif isa(self.source,"IMNeutralLabelModesBasis")
                bounds = self.source.neutralIntervals;
            else
                bounds = self.zDomain;
            end
            % Keep the original construction rule after selecting columns.
            % This preserves the rule's numerical meaning without claiming
            % that an arbitrary intrinsic evaluator is resolved by it.
            n = max(128,4*self.normalizationColumnCount_+32);
            offDiagonal = (1:n-1)./sqrt(4*(1:n-1).^2-1);
            [vectors,points] = eig(diag(offDiagonal,1)+diag(offDiagonal,-1),"vector");
            [points,index] = sort(points);
            unitWeights = 2*vectors(1,index).'.^2;
            z = zeros(0,1);
            weights = zeros(0,1);
            for j = 1:size(bounds,1)
                halfLength = diff(bounds(j,:))/2;
                z = [z;mean(bounds(j,:))+halfLength*points]; %#ok<AGROW>
                weights = [weights;halfLength*unitWeights]; %#ok<AGROW>
            end
            if includeExtrema && startsWith(self.normalization,"maxAbs")
                z = unique([z;bounds(:)]);
            end
        end
    end
end
