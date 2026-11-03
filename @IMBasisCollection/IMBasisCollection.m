classdef IMBasisCollection
    % Preserve continuous bases and their requested wavenumber-page mapping.
    %
    % A collection keeps the original scientific columns and normalization.
    % It performs no solve, truncation, acceptance decision, or approximate
    % sharing. `basisIndex` maps requested pages to stored basis objects;
    % `sourcePage` selects a page within a multi-wavenumber boundary basis.
    % Supported bases are value objects, so subsequent changes to the input
    % basis normalization do not change this collection.
    %
    % ```matlab
    % collection = IMBasisCollection({basis1,basis2},kappa=[k2 k1 k2],basisIndex=[2 1 2]);
    % G = collection.evaluate(z,variable="G",pages=[3 1]);
    % ```
    %
    % - Topic: Create collections
    % - Topic: Inspect collections
    % - Topic: Evaluate collections
    % - Declaration: classdef IMBasisCollection

    properties (SetAccess = private)
        % Continuous value bases, in stored-basis order.
        % - Topic: Inspect collections
        bases (1,:) cell
        % Requested horizontal wavenumbers, in requested-page order.
        %
        % Empty unless every page has an explicit or solved scalar identity.
        % Wavenumber-independent bases accept caller-supplied page identities;
        % per-column fixed-frequency roots cannot be replaced by one kappa.
        % Values are never rounded or replaced by representative wavenumbers.
        % - Topic: Inspect collections
        kappa (1,:) double
        % Requested-page to stored-basis mapping.
        % - Topic: Inspect collections
        basisIndex (1,:) double
        % Page within each mapped basis, normally one for scalar solves.
        % - Topic: Inspect collections
        sourcePage (1,:) double
        % Scientific column identity and evaluation capabilities per basis.
        %
        % `capabilities` preserves the authoritative common-basis record.
        % The top-level `variables` and `derivativeOrders` are evaluation
        % summaries copied from that record. Conditional physical projection
        % and sampled recovery are distinct operations; observability requires
        % preparing the actual observation system.
        % Mode labels, endpoint labels, and wavenumber pages remain separate.
        % Frequency signs are not introduced by a vertical basis collection.
        % - Topic: Inspect collections
        metadata (1,:) struct
    end

    methods (Access = private)
        [projection,recipe] = prepareProjection(self,z,weights,options)
        columns = allColumnIndices(self,page)
    end

    methods
        function self = IMBasisCollection(bases,options)
            % Construct a collection without solving or resampling modes.
            %
            % Multi-page boundary bases require explicit `sourcePage` when
            % mapping multiple requested pages. When all mapped bases have
            % known wavenumbers, omitted `kappa` is inferred from them.
            %
            % - Topic: Create collections
            % - Parameter bases: nonempty cell row of continuous value bases
            % - Parameter options.kappa: requested nonnegative wavenumber row
            % - Parameter options.basisIndex: mapping to stored bases
            % - Parameter options.sourcePage: page within each mapped basis
            % - Returns self: immutable collection of continuous bases
            arguments (Input)
                bases (1,:) cell
                options.kappa (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
                options.basisIndex (1,:) double {mustBeInteger,mustBePositive} = []
                options.sourcePage (1,:) double {mustBeInteger,mustBePositive} = []
            end
            if isempty(bases)
                error("IMBasisCollection:EmptyCollection","Supply at least one continuous basis.");
            end
            for iBasis = 1:numel(bases)
                basis = bases{iBasis};
                supported = isa(basis,"IMModeBasis") || isa(basis,"IMBasisSet") || isa(basis,"IMAnalyticalInternalModesBasis") || isa(basis,"IMGeostrophicZeroAPVModesBasis") || isa(basis,"IMAnalyticalGeostrophicZeroAPVModesBasis") || isa(basis,"IMNeutralLabelModesBasis");
                if ~supported || ~isscalar(basis) || isa(basis,"handle")
                    error("IMBasisCollection:UnsupportedBasis","Entry %d must be a supported scalar continuous value basis.",iBasis);
                end
                entry = IMBasisCollection.describeBasis(basis);
                if iBasis == 1
                    self.metadata = entry;
                else
                    self.metadata(iBasis) = entry;
                end
            end
            self.bases = bases;
            self.basisIndex = options.basisIndex;
            if isempty(self.basisIndex)
                self.basisIndex = 1:numel(bases);
            end
            if any(self.basisIndex > numel(bases))
                error("IMBasisCollection:InvalidBasisIndex","basisIndex must select existing entries in bases.");
            end
            nPages = numel(self.basisIndex);
            self.sourcePage = options.sourcePage;
            if isempty(self.sourcePage)
                self.sourcePage = ones(1,nPages);
            elseif numel(self.sourcePage) ~= nPages
                error("IMBasisCollection:InvalidSourcePage","sourcePage must contain one index per requested page.");
            end
            knownKappa = nan(1,nPages);
            for iPage = 1:nPages
                basis = bases{self.basisIndex(iPage)};
                if isa(basis,"IMModeBasis") && basis.horizontalKind == "perColumn" && ~isempty(options.kappa)
                    error("IMBasisCollection:UnsupportedHorizontalPage","A per-column dispersion basis cannot be assigned one scalar kappa. Preserve its signed kSquared roots and use a fixed-wavenumber family for scalar-wavenumber pages.");
                end
                if isa(basis,"IMModeBasis") && basis.horizontalKind == "unspecified" && ~isempty(options.kappa)
                    error("IMBasisCollection:UnsupportedHorizontalPage","A common basis with unspecified horizontal identity cannot be assigned kappa. Declare a fixed wavenumber or explicitly wavenumber-independent structure when constructing the basis.");
                end
                basisKappa = IMBasisCollection.basisWavenumbers(basis);
                if self.sourcePage(iPage) > max(1,numel(basisKappa))
                    error("IMBasisCollection:InvalidSourcePage","Requested page %d selects a nonexistent source page.",iPage);
                end
                if ~isempty(basisKappa)
                    knownKappa(iPage) = basisKappa(self.sourcePage(iPage));
                end
            end
            self.kappa = options.kappa;
            if isempty(self.kappa)
                if all(isfinite(knownKappa))
                    self.kappa = knownKappa;
                end
            elseif numel(self.kappa) ~= nPages
                error("IMBasisCollection:InvalidKappaCount","kappa must contain one value per requested page.");
            elseif any(self.kappa(isfinite(knownKappa)) ~= knownKappa(isfinite(knownKappa)))
                error("IMBasisCollection:KappaMismatch","Requested kappa must equal the underlying solved wavenumber. Approximate sharing is not supported.");
            end
        end

        values = evaluate(self,z,options)
        evaluateChunks(self,z,consumer,options)
    end

    methods (Static, Access = private)
        function entry = describeBasis(basis)
            entry = struct("family",string(class(basis)),"columnKind","mode","columnLabels",strings(1,0),"modeNumber",[],"endpointLabels",strings(1,0),"normalization","","zDomain",basis.zDomain,"variables","u","derivativeOrders",{{[0 1]}},"capabilities",struct(),"horizontalIdentity",struct(),"scientificSpecification",struct(),"constructionRecord",struct(),"rotationName","","rotationMatrix",[],"sourceCoordinateMap",[],"rotationMetadata",struct());
            if isa(basis,"IMModeBasis")
                entry.family = basis.family;
                entry.columnLabels = basis.columnLabels;
                entry.normalization = basis.normalization;
                entry.capabilities = basis.capabilities;
                entry.horizontalIdentity = struct(kind=basis.horizontalKind,k=basis.k,kSquared=basis.kSquared);
                entry.scientificSpecification = basis.scientificSpecification;
                entry.constructionRecord = basis.constructionRecord;
                entry.variables = entry.capabilities.variables;
                entry.derivativeOrders = entry.capabilities.derivativeOrders;
                if basis.family == "geostrophicZeroAPV"
                    entry.columnKind = "endpoint";
                    entry.endpointLabels = basis.canonicalEndpointLabels;
                    entry.sourceCoordinateMap = basis.sourceCoordinateMap;
                    entry.rotationMetadata = basis.rotationMetadata;
                    entry.rotationName = basis.rotationMetadata.name;
                end
                if numel(entry.columnLabels) ~= basis.nColumns
                    error("IMBasisCollection:InvalidColumnLabels","Mode labels must match the facade column count.");
                end
            elseif isa(basis,"IMNeutralLabelModesBasis")
                entry.family = "neutralLabelModes";
                entry.columnKind = "mode";
                entry.columnLabels = basis.columnLabels;
                entry.normalization = basis.normalization;
                entry.variables = ["G" "APV"];
                entry.derivativeOrders = {0:2,0};
            elseif isa(basis,"IMGeostrophicZeroAPVModesBasis") || isa(basis,"IMAnalyticalGeostrophicZeroAPVModesBasis")
                entry.family = "geostrophicZeroAPVModes";
                entry.columnKind = "endpoint";
                entry.columnLabels = reshape(string(basis.endpoints),1,[]);
                entry.endpointLabels = entry.columnLabels;
                entry.rotationName = string(basis.rotationName);
                entry.rotationMatrix = basis.rotationMatrix;
                entry.sourceCoordinateMap = basis.rotationMatrix;
                if basis.normalizationConvention ~= "unitEndpointResponse"
                    entry.columnLabels = entry.rotationName + string(1:numel(entry.endpointLabels));
                end
                entry.normalization = string(basis.normalizationConvention);
                entry.variables = ["F" "G"];
                entry.derivativeOrders = {[0 1],0};
            else
                entry.family = string(basis.evp.name);
                entry.modeNumber = reshape(basis.modeNumber,1,[]);
                entry.columnLabels = string(entry.modeNumber);
                entry.normalization = string(basis.normalization);
                if numel(entry.modeNumber) ~= numel(basis.eigenvalues)
                    error("IMBasisCollection:InvalidColumnLabels","Mode labels must match eigenvalue columns.");
                end
                if isa(basis,"IMInternalModesBasis") || isa(basis,"IMAnalyticalInternalModesBasis")
                    entry.variables = ["F" "G"];
                    entry.derivativeOrders = {0,0};
                    entry.derivativeOrders{entry.variables == basis.evp.formulation} = [0 1];
                end
            end
            if ~isa(basis,"IMModeBasis")
                entry.capabilities = IMBasisCollection.describeLegacyCapabilities(entry,basis);
                k = IMBasisCollection.basisWavenumbers(basis);
                kind = "unspecified";
                if isscalar(k), kind = "fixedWavenumber"; elseif ~isempty(k), kind = "fixedWavenumberPages"; end
                entry.horizontalIdentity = struct(kind=kind,k=k,kSquared=[]);
                entry.scientificSpecification = struct(availability="notCommonBasis");
                entry.constructionRecord = struct(availability="notCommonBasis");
            end
            supportsEmptyColumns = isa(basis,"IMModeBasis") || isa(basis,"IMNeutralLabelModesBasis");
            if (isempty(entry.columnLabels) && ~supportsEmptyColumns) || numel(unique(entry.columnLabels)) ~= numel(entry.columnLabels)
                error("IMBasisCollection:InvalidColumnLabels","Each retained column must have a distinct scientific label; only common and neutral-label bases support an empty inventory.");
            end
        end

        function description = describeLegacyCapabilities(entry,basis)
            % Keep the direct low-level collection adapters explicit while
            % the common facade supplies its own authoritative description.
            aliases = struct("variable",{},"baseVariable",{},"derivativeOffset",{});
            if ismember("F",entry.variables)
                aliases = struct(variable="Fz",baseVariable="F",derivativeOffset=1);
            end
            if ismember("G",entry.variables)
                aliases = [aliases,struct("variable",{"Gz","Gzz"},"baseVariable",{"G","G"},"derivativeOffset",{1,2})];
            end
            spectra = strings(1,0);
            for name = ["h","eigenvalues"]
                if isprop(basis,name), spectra(end+1) = name; end %#ok<AGROW>
            end
            normalizations = entry.normalization;
            if ismethod(basis,"normalizationNames"), normalizations = basis.normalizationNames(); end
            if isa(basis,"IMNeutralLabelModesBasis")
                derivativeConvention = "classicalIntervalInterior";
                traceConvention = "insideIntervalEndpoints";
                operations = struct(variable="APV",kind="prescribedContinuousDual",availability="conditional",requirements=["APV observations with fixed dz quadrature","nonzero f0 and representable continuous APV norms for retained columns","empty inventories preserve an empty physical dual"]);
            else
                derivativeConvention = "classicalPhysical";
                traceConvention = "physicalDomainEvaluation";
                operations = struct("variable",{},"kind",{},"availability",{},"requirements",{});
                for variable = entry.variables
                    operations(end+1) = struct(variable=variable,kind="familyRecipe",availability="conditional",requirements="The requested variable must have an available family-owned physical pairing recipe."); %#ok<AGROW>
                end
            end
            recovery = struct(synthesisVariables=entry.variables,preparation="IMDiscreteTransform.fromObservationMatrices",availability="requiresPreparedSystem",requirements=["explicit observation model and column identities","explicit scientific column activity","declared fitting or observation factor, or physical dual","successful rank and overlap assessment of the prepared system"],noiseEvidence="requiresDeclaredObservationFactor");
            description = struct(variables=entry.variables,derivativeOrders={entry.derivativeOrders},aliases=aliases,derivativeConvention=derivativeConvention,traceConvention=traceConvention,traceOrders={entry.derivativeOrders},spectrumNames=spectra,normalizationNames=normalizations,normalizationRequirements="Each normalization rule validates its retained-column magnitudes.",projectionOperations=operations,sampledRecovery=recovery,supportsEmptyColumns=isa(basis,"IMNeutralLabelModesBasis"));
        end

        function kappa = basisWavenumbers(basis)
            if isa(basis,"IMModeBasis")
                kappa = [];
                if basis.horizontalKind == "fixedWavenumber", kappa = basis.k; end
            elseif isprop(basis,"k")
                kappa = reshape(basis.k,1,[]);
            elseif isprop(basis,"evp") && isfield(basis.evp.parameters,"k")
                kappa = reshape(basis.evp.parameters.k,1,[]);
            else
                kappa = [];
            end
        end
    end
end
