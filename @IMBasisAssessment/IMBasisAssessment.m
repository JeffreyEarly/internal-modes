classdef IMBasisAssessment
    % Measure a fixed projection without selecting or replacing its columns.
    %
    % Gram, leakage, and supplied-product errors share one result vocabulary.
    % `measurements` contains one row per requested prefix and quantity.
    % Measurement status is measured, inconclusive, or notRequested; an
    % unqualified product reference is inconclusive even if its error is small.
    % Apply explicit tolerances with `applyPolicy` to obtain a separate decision.
    %
    % Product recipes supply signed continuous pairings, not full-band
    % coefficients: each prefix solves its own continuous target Gram system.
    % The sampled projection keeps its signed metric. Positive coefficient
    % error norms use the projection's separately supplied majorant.
    %
    % ```matlab
    % assessment = IMBasisAssessment(projection,prefixCounts=1:8);
    % decision = assessment.applyPolicy(gramTolerance=1e-2);
    % ```
    %
    % - Topic: Measure projections
    % - Topic: Apply acceptance policies
    % - Declaration: classdef IMBasisAssessment

    properties (SetAccess = private)
        % Original fixed projection, with every requested column preserved.
        % - Topic: Measure projections
        projection
        % Scientific family, variable, page, and other caller-supplied identity.
        % - Topic: Measure projections
        identity
        % Mode columns or endpoint coordinates; only modes have prefixes.
        % - Topic: Measure projections
        columnKind
        % Explicitly examined column counts, including the full requested band.
        % - Topic: Measure projections
        prefixCounts
        % Quantity, value, status, limiting inputs, and reference provenance.
        % - Topic: Measure projections
        measurements
        % Actual sampled coverage and omissions, without exhaustive guarantees.
        % - Topic: Measure projections
        coverage
        % Measured assessment time, separate from caller construction costs.
        % - Topic: Measure projections
        costs
    end

    methods
        function self = IMBasisAssessment(projection,options)
            % Measure Gram quality and optional prepared leakage and products.
            %
            % `leakage` contains sampledValues, normSquared, and columnLabels.
            % Check columns must belong to this same scientific family and page.
            % Labels matching retained columns are excluded at each prefix.
            % `products` contains sampledValues, referencePairings,
            % normSquared, inputColumns (2-by-nProducts in the selected output
            % column order; zero denotes an input
            % outside this output family), inputLabels (2-by-nProducts), and
            % referenceStatus (qualified or inconclusive/unverified).
            % Optional referenceProvenance describes the independent evidence.
            %
            % - Topic: Measure projections
            % - Declaration: assessment = IMBasisAssessment(projection,options)
            arguments
                projection (1,1) IMProjection
                options.identity (1,1) struct = struct()
                options.columnKind (1,1) string {mustBeMember(options.columnKind,["mode","endpoint"])} = "mode"
                options.prefixCounts (1,:) double {mustBeInteger,mustBePositive} = []
                options.leakage struct = struct.empty
                options.products struct = struct.empty
                options.coverage (1,1) struct = struct()
                options.constructionSeconds (1,1) double = NaN
            end
            started = tic;
            nColumns = projection.columnCount;
            if nColumns == 0
                error("IMBasisAssessment:EmptyProjection","Supply a projection with at least one scientific column.");
            end
            counts = options.prefixCounts;
            if isempty(counts)
                counts = nColumns;
            end
            if any(counts > nColumns) || numel(unique(counts)) ~= numel(counts) || any(diff(counts) <= 0) || counts(end) ~= nColumns
                error("IMBasisAssessment:InvalidPrefixes","prefixCounts must increase without duplicates, end at the requested column count, and stay within that count.");
            end
            if options.columnKind == "endpoint" && ~isequal(counts,nColumns)
                error("IMBasisAssessment:EndpointPrefix","Endpoint coordinates are assessed as the complete requested set, never as modal prefixes.");
            end
            leakage = options.leakage;
            products = options.products;
            if ~isempty(leakage)
                requireFields(leakage,["sampledValues","normSquared","columnLabels"],"leakage");
                if size(leakage.sampledValues,1) ~= projection.sampleCount || numel(leakage.normSquared) ~= size(leakage.sampledValues,2) || numel(leakage.columnLabels) ~= size(leakage.sampledValues,2)
                    error("IMBasisAssessment:InvalidLeakage","Leakage samples, norms, and labels must describe the same check columns.");
                end
                leakage.columnLabels = reshape(string(leakage.columnLabels),1,[]);
                leakage.normSquared = reshape(leakage.normSquared,1,[]);
                if any(~isfinite(leakage.normSquared) | leakage.normSquared <= 0) || numel(unique(leakage.columnLabels)) ~= numel(leakage.columnLabels)
                    error("IMBasisAssessment:InvalidLeakage","Leakage check columns require distinct labels and finite positive norm squares.");
                end
            end
            if ~isempty(products)
                requireFields(products,["sampledValues","referencePairings","normSquared","inputColumns","inputLabels","referenceStatus"],"products");
                nProducts = size(products.sampledValues,2);
                if size(products.sampledValues,1) ~= projection.sampleCount || ~isequal(size(products.referencePairings),[nColumns nProducts]) || numel(products.normSquared) ~= nProducts || ~isequal(size(products.inputColumns),[2 nProducts]) || ~isequal(size(products.inputLabels),[2 nProducts])
                    error("IMBasisAssessment:InvalidProducts","Product samples, reference pairings, norms, and input identities must describe the same products and requested output columns.");
                end
                if any(~isfinite(products.inputColumns) | products.inputColumns < 0 | products.inputColumns > nColumns | fix(products.inputColumns) ~= products.inputColumns,"all")
                    error("IMBasisAssessment:InvalidProducts","Product inputColumns must contain ordinals within the output family; zero denotes an input outside this output family.");
                end
                if any(~isfinite(products.sampledValues),"all") || any(~isfinite(products.referencePairings),"all") || ~isreal(products.normSquared) || any(~isfinite(products.normSquared) | products.normSquared < 0,"all")
                    error("IMBasisAssessment:InvalidProducts","Products and signed reference pairings must be finite; norm squares must be real, finite, and nonnegative.");
                end
                products.normSquared = reshape(products.normSquared,1,[]);
                products.inputLabels = string(products.inputLabels);
                products.referenceStatus = string(products.referenceStatus);
                if ~isscalar(products.referenceStatus) || ~ismember(products.referenceStatus,["qualified","inconclusive","unverified"])
                    error("IMBasisAssessment:InvalidReferenceStatus","Product referenceStatus must be qualified, inconclusive, or unverified.");
                end
                if ~isfield(products,"referenceProvenance")
                    products.referenceProvenance = "unspecified";
                end
            end
            nRows = 3*numel(counts);
            columnCount = reshape(repelem(counts(:),3),[],1);
            quantity = repmat(["gram";"leakage";"quadraticAliasing"],numel(counts),1);
            value = nan(nRows,1);
            status = repmat("notRequested",nRows,1);
            limitingInputI = strings(nRows,1);
            limitingInputJ = strings(nRows,1);
            referenceStatus = repmat("notRequired",nRows,1);
            referenceProvenance = strings(nRows,1);
            examinedCount = zeros(nRows,1);
            for iPrefix = 1:numel(counts)
                count = counts(iPrefix);
                p = projection.prefix(count);
                row = 3*(iPrefix-1)+1;
                value(row) = p.gramError;
                status(row) = "measured";
                if ~p.supportsGramAssessment
                    status(row) = "unsupported";
                end
                examinedCount(row) = count;
                if ~isempty(leakage)
                    rejected = ~ismember(leakage.columnLabels,p.columnLabels);
                    examinedCount(row+1) = nnz(rejected);
                    if any(rejected)
                        errors = p.leakage(leakage.sampledValues(:,rejected),leakage.normSquared(rejected));
                        [value(row+1),limiting] = max(errors);
                        labels = leakage.columnLabels(rejected);
                        limitingInputI(row+1) = labels(limiting);
                        status(row+1) = "measured";
                    else
                        status(row+1) = "inconclusive";
                    end
                end
                if ~isempty(products)
                    selected = all(products.inputColumns <= count,1);
                    examinedCount(row+2) = nnz(selected);
                    referenceStatus(row+2) = products.referenceStatus;
                    referenceProvenance(row+2) = string(products.referenceProvenance);
                    if any(selected)
                        active = p.activeColumnMask;
                        coefficients = zeros(count,nnz(selected));
                        pairings = products.referencePairings(1:count,selected);
                        if any(pairings(~active,:) ~= 0,"all")
                            error("IMBasisAssessment:InactiveReferencePairing","Inactive output columns require zero reference pairings; supplied reference data are inconsistent with the projection.");
                        end
                        target = p.targetGramMatrix(active,active);
                        if rank(target) < nnz(active)
                            error("IMBasisAssessment:SingularReferenceTarget","Continuous reference projection requires a nonsingular target on active columns.");
                        end
                        coefficients(active,:) = target \ pairings(active,:);
                        errors = p.productError(products.sampledValues(:,selected),coefficients,products.normSquared(selected));
                        [value(row+2),limiting] = max(errors);
                        labels = products.inputLabels(:,selected);
                        limitingInputI(row+2) = labels(1,limiting);
                        limitingInputJ(row+2) = labels(2,limiting);
                        status(row+2) = "measured";
                    else
                        status(row+2) = "inconclusive";
                    end
                    if products.referenceStatus ~= "qualified"
                        status(row+2) = "inconclusive";
                    end
                end
            end
            self.projection = projection;
            self.identity = options.identity;
            self.columnKind = options.columnKind;
            self.prefixCounts = counts;
            self.measurements = table(columnCount,quantity,value,status,limitingInputI,limitingInputJ,referenceStatus,referenceProvenance,examinedCount);
            self.coverage = options.coverage;
            self.coverage.exhaustiveGuarantee = false;
            self.coverage.superpositionGuarantee = false;
            self.costs = struct(constructionSeconds=options.constructionSeconds,assessmentSeconds=toc(started));
        end

        function decision = applyPolicy(self,options)
            % Apply explicit tolerances without changing the assessed basis.
            %
            % Disabled quantities do not enter the decision. Missing or
            % unqualified enabled measurements are inconclusive. A measured
            % failure rejects even when another enabled quantity is inconclusive.
            % Prefix recommendations are diagnostic only; requestedColumnCount
            % and the stored projection are never reduced.
            %
            % - Topic: Apply acceptance policies
            % - Declaration: decision = applyPolicy(assessment,options)
            arguments
                self (1,1) IMBasisAssessment
                options.gramTolerance double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
                options.leakageTolerance double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
                options.quadraticAliasingTolerance double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
            end
            names = ["gram","leakage","quadraticAliasing"];
            tolerances = {options.gramTolerance,options.leakageTolerance,options.quadraticAliasingTolerance};
            if any(cellfun(@numel,tolerances) > 1) || all(cellfun(@isempty,tolerances))
                error("IMBasisAssessment:InvalidPolicy","Supply at least one scalar tolerance; omit quantities that are not acceptance requirements.");
            end
            outcomes = repmat("accepted",numel(self.prefixCounts),1);
            for iPrefix = 1:numel(self.prefixCounts)
                for iQuantity = 1:numel(names)
                    if isempty(tolerances{iQuantity})
                        continue
                    end
                    row = self.measurements(self.measurements.columnCount == self.prefixCounts(iPrefix) & self.measurements.quantity == names(iQuantity),:);
                    if row.status ~= "measured" || isnan(row.value)
                        if outcomes(iPrefix) ~= "rejected"
                            outcomes(iPrefix) = "inconclusive";
                        end
                    elseif row.value > tolerances{iQuantity}
                        outcomes(iPrefix) = "rejected";
                    end
                end
            end
            consecutive = cumprod(double(outcomes == "accepted")) ~= 0;
            acceptedCounts = self.prefixCounts(consecutive);
            if isempty(acceptedCounts)
                largestExaminedAcceptedPrefix = 0;
            else
                largestExaminedAcceptedPrefix = acceptedCounts(end);
            end
            overallStatus = outcomes(end);
            if any(outcomes == "rejected")
                overallStatus = "rejected";
            elseif any(outcomes == "inconclusive")
                overallStatus = "inconclusive";
            end
            decision = struct(status=overallStatus,requestedBandStatus=outcomes(end),requestedColumnCount=self.projection.columnCount,largestExaminedAcceptedPrefix=largestExaminedAcceptedPrefix,allPrefixesExamined=isequal(self.prefixCounts,1:self.projection.columnCount),prefixDecisions=table(self.prefixCounts(:),outcomes,consecutive,VariableNames=["columnCount","status","consecutiveAccepted"]),tolerances=options);
        end
    end
end

function requireFields(value,fields,name)
if ~isscalar(value) || ~all(isfield(value,fields))
    error("IMBasisAssessment:InvalidPreparedData","%s must be a scalar struct containing %s.",name,join(fields,", "));
end
end
