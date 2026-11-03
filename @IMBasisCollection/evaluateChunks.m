function evaluateChunks(self,z,consumer,options)
% Stream continuous basis values in bounded sample and page chunks.
%
% Calls `consumer(values,rowIndices,pagePositions,requestedPages)`, where
% pagePositions indexes the supplied `pages` vector and requestedPages holds
% the corresponding collection page IDs. This preserves repeated page IDs
% without ambiguity. Values have shape nRows x nColumns x nChunkPages.
% No full output array is allocated. Validation precedes the first callback.
% Compatible native spectral bases share evaluation matrices exactly;
% analytical bases and custom subclasses retain their own evaluators.
% Custom evaluators control their internal allocations; chunk bounds apply
% to callbacks and the built-in evaluation paths, not arbitrary user code.
% Neutral-label G derivatives are interval-interior derivatives with inside
% endpoint traces. In particular, order 2 omits distributional interface terms.
% A normalization override applies to every selected basis for this call.
% Omit it to preserve each basis's stored normalization. Fixed endpoint
% conventions accept only their existing normalization; use an explicit
% boundary-coordinate rotation to change their scientific columns.
%
% - Topic: Evaluate collections
% - Parameter z: physical coordinate column
% - Parameter consumer: callback receiving values, rows, positions, and page IDs
% - Parameter options.variable: required scalar variable name
% - Parameter options.derivativeOrder: supported vertical derivative order
% - Parameter options.normalization: optional common normalization override
% - Parameter options.pages: requested collection pages, empty means all
% - Parameter options.columns: explicit scientific column indices
% - Parameter options.sampleChunkSize: maximum coordinate rows per callback
% - Parameter options.pageChunkSize: maximum requested pages per callback
arguments (Input)
    self (1,1) IMBasisCollection
    z (:,1) double {mustBeReal,mustBeFinite}
    consumer (1,1) function_handle
    options.variable (1,1) string = ""
    options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.normalization = []
    options.pages (1,:) double {mustBeInteger,mustBePositive} = []
    options.columns (1,:) double {mustBeInteger,mustBePositive} = []
    options.sampleChunkSize (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 1024
    options.pageChunkSize (1,1) double {mustBeInteger,mustBePositive,mustBeFinite} = 16
end
pages = options.pages;
if isempty(pages)
    pages = 1:numel(self.basisIndex);
elseif any(pages > numel(self.basisIndex))
    error("IMBasisCollection:InvalidPage","pages must select existing requested pages.");
end
variable = options.variable;
if variable == ""
    error("IMBasisCollection:MissingVariable","Supply variable explicitly, such as variable=""G"". Inspect metadata.variables for supported channels.");
end
columnCounts = arrayfun(@(i) numel(self.metadata(self.basisIndex(i)).columnLabels),pages);
columns = options.columns;
if isempty(columns)
    if any(columnCounts ~= columnCounts(1))
        error("IMBasisCollection:HeterogeneousColumns","Selected pages have different column counts. Evaluate homogeneous page groups or select explicit columns.");
    end
    columns = 1:columnCounts(1);
elseif any(columns > min(columnCounts))
    error("IMBasisCollection:InvalidColumn","columns must exist on every selected page. Evaluate differing column counts separately.");
end
indices = unique(self.basisIndex(pages),"stable");
factors = cell(size(self.bases));
for iBasis = indices
    entry = self.metadata(iBasis);
    iVariable = find(entry.variables == variable,1);
    if isempty(iVariable)
        error("IMBasisCollection:UnsupportedVariable","Basis %d does not support %s. Inspect metadata.variables.",iBasis,variable);
    end
    if ~ismember(options.derivativeOrder,entry.derivativeOrders{iVariable})
        error("IMBasisCollection:UnsupportedDerivative","Basis %d does not implement derivative order %d for %s.",iBasis,options.derivativeOrder,variable);
    end
    if any(z < entry.zDomain(1) | z > entry.zDomain(2))
        error("IMBasisCollection:OutsideDomain","Evaluation coordinates must lie inside every selected basis domain.");
    end
    basis = self.bases{iBasis};
    if entry.columnKind == "endpoint"
        if ~isempty(options.normalization) && (~isscalar(string(options.normalization)) || string(options.normalization) ~= string(basis.normalizationConvention))
            error("IMBasisCollection:UnsupportedNormalization","Basis %d has fixed normalization %s. Change endpoint coordinates with an explicit basis rotation before constructing the collection.",iBasis,basis.normalizationConvention);
        end
    else
        if isempty(options.normalization)
            normalization = basis.normalizationFactors();
        else
            normalization = basis.normalizationFactors(options.normalization);
        end
        factors{iBasis} = normalization(columns);
    end
end
% Normalize each immutable stored basis once, including repeated requests.
[groups,basisGroup] = prepareNativeGroups(self,indices,variable,options.derivativeOrder);
normalizationArguments = {};
if ~isempty(options.normalization)
    normalizationArguments = namedargs2cell(struct(normalization=options.normalization));
end
% Sample chunks are outermost so evaluation matrices can serve every page
% chunk without retaining operators for the complete physical grid.
for firstRow = 1:options.sampleChunkSize:numel(z)
    rows = firstRow:min(numel(z),firstRow+options.sampleChunkSize-1);
    zChunk = z(rows);
    matrices = cell(size(groups));
    for firstPage = 1:options.pageChunkSize:numel(pages)
        positions = firstPage:min(numel(pages),firstPage+options.pageChunkSize-1);
        chunkPages = pages(positions);
        storedIndices = unique(self.basisIndex(chunkPages),"stable");
        values = zeros(numel(rows),numel(columns),numel(positions));
        for iGroup = unique(basisGroup(storedIndices))
            if iGroup == 0, continue; end
            group = groups{iGroup};
            members = storedIndices(basisGroup(storedIndices)==iGroup);
            if isempty(matrices{iGroup})
                identity = eye(group.nPolynomials);
                if group.derivativeOrder == 0
                    matrices{iGroup} = group.solver.evaluateNativeModes(identity,zChunk);
                else
                    matrices{iGroup} = group.solver.evaluatePhysicalDerivative(identity,zChunk,group.derivativeOrder);
                end
            end
            native = cell(1,numel(members));
            for iMember = 1:numel(members)
                native{iMember} = self.bases{members(iMember)}.nativeModes(:,columns);
            end
            sampled = matrices{iGroup}*cat(2,native{:});
            for iMember = 1:numel(members)
                iBasis = members(iMember);
                memberColumns = (iMember-1)*numel(columns)+(1:numel(columns));
                memberValues = sampled(:,memberColumns);
                basis = self.bases{iBasis};
                if options.derivativeOrder == 0 && isa(basis,"IMInternalModesBasis") && variable ~= string(basis.evp.formulation)
                    context = basis.evp.contextForSolver(basis.solver);
                    if variable == "F"
                        memberValues = basis.evp.FfromGz(zChunk,memberValues,basis.h(columns),context);
                    else
                        memberValues = basis.evp.GfromFz(zChunk,memberValues,basis.h(columns),context);
                    end
                end
                memberValues = memberValues./factors{iBasis};
                selected = find(self.basisIndex(chunkPages) == iBasis);
                values(:,:,selected) = repmat(memberValues,1,1,numel(selected));
            end
        end
        for iBasis = storedIndices(basisGroup(storedIndices)==0)
            basis = self.bases{iBasis};
            entry = self.metadata(iBasis);
            selected = find(self.basisIndex(chunkPages) == iBasis);
            if entry.columnKind == "endpoint"
                sourcePages = self.sourcePage(chunkPages(selected));
                sampled = basis.evaluate(zChunk,normalizationArguments{:},variable=variable,derivativeOrder=options.derivativeOrder,pages=sourcePages);
                sampled = sampled(:,columns,:);
            else
                sampled = basis.evaluate(zChunk,normalizationArguments{:},variable=variable,derivativeOrder=options.derivativeOrder);
                if size(sampled,2) ~= numel(entry.columnLabels)
                    error("IMBasisCollection:InvalidEvaluationShape","Basis evaluation must return one column per scientific label.");
                end
                sampled = repmat(sampled(:,columns),1,1,numel(selected));
            end
            if size(sampled,1) ~= numel(rows) || size(sampled,2) ~= numel(columns) || size(sampled,3) ~= numel(selected)
                error("IMBasisCollection:InvalidEvaluationShape","Basis evaluation must match the requested row, column, and page counts.");
            end
            values(:,:,selected) = sampled;
        end
        consumer(values,rows,positions,chunkPages);
    end
end
end

function [groups,basisGroup] = prepareNativeGroups(self,indices,variable,derivativeOrder)
groups = {};
basisGroup = zeros(1,numel(self.bases));
for iBasis = indices
    basis = self.bases{iBasis};
    if ~ismember(string(class(basis)),["IMBasisSet","IMInternalModesBasis"]) || ~strcmp(class(basis.solver),'IMSolverSpectral')
        continue
    end
    nativeOrder = derivativeOrder;
    if isa(basis,"IMInternalModesBasis") && variable ~= string(basis.evp.formulation)
        if derivativeOrder > 0
            continue
        end
        nativeOrder = 1;
    end
    nPolynomials = size(basis.nativeModes,1);
    iGroup = find(cellfun(@(g) g.derivativeOrder == nativeOrder && g.nPolynomials == nPolynomials && isequal(g.solver,basis.solver),groups),1);
    if isempty(iGroup)
        groups{end+1} = struct("solver",basis.solver,"derivativeOrder",nativeOrder,"nPolynomials",nPolynomials); %#ok<AGROW>
        iGroup = numel(groups);
    end
    basisGroup(iBasis) = iGroup;
end
end
