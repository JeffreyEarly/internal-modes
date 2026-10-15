classdef IMModeBasis
    % Evaluate one mode family with shared identities and normalization names.
    %
    % Factories in `IMModes` create these immutable continuous bases. Wave
    % profiles obey $$F_j = D G_j'$$, where $$D = z_s-z_b$$. Every
    % normalization multiplies all fields in one column by the same
    % positive factor; signed metric signatures remain explicit.
    %
    % ```matlab
    % basis = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=8);
    % F = basis.F(z);
    % ```
    %
    % - Topic: Inspect modes
    % - Topic: Evaluate modes
    % - Topic: Project modes
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
        % Available profile evaluation names.
        % - Topic: Inspect modes
        capabilities (1,:) string
        % Immutable physical column configuration.
        % - Topic: Inspect modes
        background IMModeBackground
        % Increasing bottom and surface coordinates, in meters.
        % - Topic: Inspect modes
        zDomain (1,2) double
        % Signs remaining after signed normalization, otherwise positive one.
        % - Topic: Inspect modes
        normSignatures (1,:) double
        % Horizontal wavenumber used in fixed-wavenumber wave normalization.
        % - Topic: Inspect modes
        k (1,1) double
        % Physical modal frequencies when defined by the family.
        % - Topic: Inspect modes
        frequency (1,:) double
        % Additional family and construction information.
        % - Topic: Inspect modes
        metadata (1,1) struct
    end

    properties (Dependent, SetAccess = private)
        % Equivalent depths; unavailable for families without equivalent depth.
        % - Topic: Inspect modes
        h
        % EVP eigenvalues; unavailable for families without an EVP.
        % - Topic: Inspect modes
        eigenvalues
    end

    properties (Access = private)
        source
        columns (1,:) double
        waveConvention (1,1) logical
        intrinsicEvaluator
        equivalentDepths_ (1,:) double
        eigenvalues_ (1,:) double
        hasEquivalentDepths_ (1,1) logical = false
        hasEigenvalues_ (1,1) logical = false
        % Immutable derived scales depend only on source, columns and normalization.
        scaleFactors_ (1,:) double
        sourceUnityFactors_ (1,:) double
        integrationIntervals_ (:,2) double
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
            % - Parameter options.source: existing numerical or analytical basis, or empty
            % - Parameter options.background: common physical column configuration
            % - Parameter options.family: scientific family name
            % - Parameter options.columns: selected source columns; all by default
            % - Parameter options.modeInfo: table with columnLabel, branch and branchIndex
            % - Parameter options.normalization: resolved common normalization name
            % - Parameter options.waveConvention: whether to enforce F = D G prime
            % - Parameter options.k: nonnegative fixed horizontal wavenumber
            % - Parameter options.h: optional equivalent depths for selected columns
            % - Parameter options.eigenvalues: optional selected eigenvalues
            % - Parameter options.hasEquivalentDepths: advertise an empty equivalent-depth spectrum
            % - Parameter options.hasEigenvalues: advertise an empty EVP spectrum
            % - Parameter options.frequency: optional physical modal frequencies
            % - Parameter options.intrinsicEvaluator: function of variable, z and derivative order
            % - Parameter options.capabilities: available variables for an intrinsic evaluator
            % - Parameter options.variables: optional intrinsic variable list
            % - Parameter options.metadata: additional scientific metadata
            % - Parameter options.integrationIntervals: exact intrinsic support intervals
            % - Returns self: immutable continuous family basis
            % - Developer: true
            arguments (Input)
                options.source = []
                options.background (1,1) IMModeBackground
                options.family (1,1) string
                options.columns (1,:) double {mustBeInteger,mustBePositive} = []
                options.modeInfo table = table()
                options.normalization (1,1) string = "depthRMSF"
                options.waveConvention (1,1) logical = false
                options.k (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 0
                options.h (1,:) double {mustBeReal} = []
                options.eigenvalues (1,:) double {mustBeReal} = []
                options.hasEquivalentDepths (1,1) logical = false
                options.hasEigenvalues (1,1) logical = false
                options.frequency (1,:) double {mustBeReal} = []
                options.intrinsicEvaluator = []
                options.capabilities (1,:) string = ["F","G","Gz","Gzz"]
                options.variables (1,:) string = strings(1,0)
                options.metadata (1,1) struct = struct()
                options.integrationIntervals (:,2) double {mustBeReal,mustBeFinite} = zeros(0,2)
            end
            arguments (Output)
                self (1,1) IMModeBasis
            end
            self.source = options.source;
            self.background = options.background;
            self.family = options.family;
            self.zDomain = options.background.zDomain;
            self.normalization = options.normalization;
            self.waveConvention = options.waveConvention;
            self.k = options.k;
            self.frequency = options.frequency;
            self.metadata = options.metadata;
            if self.waveConvention
                self.metadata.profileConvention = "F=DGz";
            end
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
                self.capabilities = options.capabilities;
                if ~isempty(options.variables)
                    self.capabilities = options.variables;
                    if ismember("G",self.capabilities)
                        self.capabilities = unique([self.capabilities,"Gz","Gzz"],"stable");
                    end
                end
                if nSource > 0 && isempty(self.intrinsicEvaluator)
                    error("IMModeBasis:MissingEvaluator","Nonempty intrinsic modes require an intrinsicEvaluator.");
                end
            elseif isa(self.source,"IMNeutralLabelModesBasis")
                nSource = numel(self.source.columnLabels);
                self.capabilities = ["G","Gz","Gzz","APV"];
            elseif isa(self.source,"IMBasisSet") || isa(self.source,"IMAnalyticalInternalModesBasis")
                nSource = numel(self.source.eigenvalues);
                self.capabilities = ["F","G","Gz","Gzz"];
                self.sourceUnityFactors_ = self.source.normalizationFactors("unity");
                factors = self.sourceUnityFactors_;
                self.source = self.source.addNormalization("unity",@(~,j) factors(j));
                self.source.normalization = "unity";
            elseif isa(self.source,"IMGeostrophicZeroAPVModesBasis") || isa(self.source,"IMAnalyticalGeostrophicZeroAPVModesBasis")
                if numel(self.source.k) ~= 1
                    error("IMModeBasis:MultipleWavenumbers","Wrap one zero-APV wavenumber page at a time.");
                end
                nSource = numel(self.source.endpoints);
                self.capabilities = ["F","G","Gz","Gzz"];
            else
                error("IMModeBasis:UnsupportedSource","source must be a supported continuous mode basis.");
            end
            if ~isempty(self.source) && ~isequal(self.source.zDomain,self.zDomain)
                error("IMModeBasis:IncompatibleDomain","The source and background must have the same zDomain.");
            end
            self.columns = options.columns;
            if isempty(self.columns)
                self.columns = 1:nSource;
            elseif any(self.columns > nSource) || numel(unique(self.columns)) ~= numel(self.columns)
                error("IMModeBasis:InvalidColumns","Select distinct existing source columns.");
            end
            self.nColumns = numel(self.columns);
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
                self.equivalentDepths_ = self.source.h(self.columns);
                self.hasEquivalentDepths_ = true;
            end
            if ~isempty(options.eigenvalues)
                self.eigenvalues_ = options.eigenvalues;
                self.hasEigenvalues_ = true;
            elseif ~isempty(self.source) && isprop(self.source,"eigenvalues")
                self.eigenvalues_ = self.source.eigenvalues(self.columns);
                self.hasEigenvalues_ = true;
            end
            if (self.hasEquivalentDepths_ && numel(self.equivalentDepths_) ~= self.nColumns) || (self.hasEigenvalues_ && numel(self.eigenvalues_) ~= self.nColumns) || (~isempty(self.frequency) && numel(self.frequency) ~= self.nColumns)
                error("IMModeBasis:InvalidSpectrumCount","h, eigenvalues and frequency must have one value per selected column when supplied.");
            end
            if self.waveConvention && ~isempty(self.source) && self.source.evp.formulation ~= "G"
                error("IMModeBasis:WaveFormulation","The standard wave convention requires a G-form source.");
            end
            [self.scaleFactors_,self.normSignatures] = self.normalizationScales();
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
            % - Parameter variable: F, G, Gz, Gzz or APV, when available
            % - Parameter derivativeOrder: physical derivative order, zero by default
            % - Returns values: normalized sample-by-column matrix
            arguments (Input)
                self (1,1) IMModeBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                variable (1,1) string {mustBeMember(variable,["F","G","Gz","Gzz","APV"])}
                derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
            end
            arguments (Output)
                values (:,:) double
            end
            if variable == "Gz"
                variable = "G";
                derivativeOrder = derivativeOrder+1;
            elseif variable == "Gzz"
                variable = "G";
                derivativeOrder = derivativeOrder+2;
            end
            if ~ismember(variable,self.capabilities)
                error("IMModeBasis:UnavailableVariable","The %s family has no %s profile.",self.family,variable);
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

        function projection = projectionOnGrid(self,z,weights,options)
            % Construct the prescribed continuous neutral-label APV dual.
            %
            % APV cannot identify a constant neutral label. Select observable
            % columns explicitly when the complete basis includes that label.
            %
            % - Topic: Project modes
            % - Parameter z: increasing physical sample column
            % - Parameter weights: fixed physical dz quadrature weights
            % - Parameter options.columns: selected facade columns, all by default
            % - Returns projection: map from APV samples to normalized-basis coefficients
            arguments (Input)
                self (1,1) IMModeBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                weights (:,1) double {mustBeReal,mustBeFinite}
                options.columns (1,:) double {mustBeInteger,mustBePositive} = 1:self.nColumns
            end
            arguments (Output)
                projection (1,1) IMProjection
            end
            if ~isa(self.source,"IMNeutralLabelModesBasis")
                error("IMModeBasis:UnavailableProjection","This basis does not provide a prescribed APV projection.");
            end
            selectedColumns = options.columns;
            if isempty(selectedColumns)
                selectedColumns = 1:self.nColumns;
            end
            if any(selectedColumns > self.nColumns) || numel(unique(selectedColumns)) ~= numel(selectedColumns)
                error("IMModeBasis:InvalidColumns","Select distinct existing basis columns.");
            end
            base = self.source.projectionOnGrid(z,weights,columns=self.columns(selectedColumns));
            factors = self.scaleFactors_(selectedColumns);
            forward = factors(:).*base.forwardMatrix;
            majorant = base.majorantGramMatrix./(factors(:)*factors);
            provenance = base.provenance;
            provenance.normalization = self.normalization;
            projection = IMProjection.fromPrescribedDual(forward,eye(numel(selectedColumns)),eye(numel(selectedColumns)),majorantGramMatrix=majorant,columnLabels=self.columnLabels(selectedColumns),provenance=provenance);
        end
    end

    methods (Hidden)
        function recipe = projectionRecipe(self,options)
            % Describe the continuous neutral-label APV dual and its reference.
            % - Topic: Project modes
            % - Parameter options.variable: APV, or empty for the family default
            % - Returns recipe: family-owned reference and capability metadata
            % - Developer: true
            arguments (Input)
                self (1,1) IMModeBasis
                options.variable (1,1) string = ""
            end
            arguments (Output)
                recipe (1,1) struct
            end
            if ~isa(self.source,"IMNeutralLabelModesBasis") || ~ismember(options.variable,["","APV"])
                error("IMModeBasis:UnavailableProjection","This basis does not provide the requested prescribed APV projection.");
            end
            recipe = self.source.projectionRecipe(variable="APV");
            factors = self.scaleFactors_;
            recipe.columnLabels = self.columnLabels;
            recipe.normalization = self.normalization;
            recipe.targetGramMatrix = eye(self.nColumns);
            recipe.majorantGramMatrix = recipe.majorantGramMatrix(self.columns,self.columns)./(factors(:)*factors);
            recipe.evaluateFunction = @(z) self.APV(z);
            recipe.provenance.normalization = self.normalization;
        end
    end

    methods (Access = private)
        function values = unscaledValues(self,variable,z,derivativeOrder)
            if self.nColumns == 0
                values = zeros(numel(z),0);
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
            values = values(:,self.columns);
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
            names = ["depthRMSF","depthRMSG","waveEnergy","generalizedEnergy","generalizedPotentialEnstrophy","endpointResponse","maxAbsF","maxAbsG"];
            if ~ismember(self.normalization,names)
                error("IMModeBasis:UnknownNormalization","Choose one of the common normalization names: %s.",join(names,", "));
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
                case "generalizedEnergy"
                    if isempty(self.sourceUnityFactors_) || ~ismethod(self.source,"gramMatrix") || self.waveConvention
                        error("IMModeBasis:UnavailableNormalization","This family has no generalized-energy G pairing.");
                    end
                    gram = self.source.gramMatrix(variable="G");
                    squares = reshape(diag(gram),1,[]);
                    squares = squares(self.columns);
                    signatures = sign(squares);
                    scales = sqrt(abs(squares));
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
                        scales(self.source.m_n(self.columns) == 0) = 1;
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
            if ~ismember(variable,self.capabilities)
                error("IMModeBasis:UnavailableNormalization","%s requires a %s profile, which this family does not provide.",self.normalization,variable);
            end
        end

        function [z,weights] = normalizationQuadrature(self)
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
                if startsWith(self.normalization,"maxAbs")
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
            n = max(128,4*self.nColumns+32);
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
            if startsWith(self.normalization,"maxAbs")
                z = unique([z;bounds(:)]);
            end
        end
    end
end
