function values = evaluate(self,z,options)
% Evaluate continuous structures on selected requested pages.
%
% The output is `nZ x nColumns x nSelectedPages` (MATLAB omits trailing
% singleton dimensions). Empty `pages` and `columns` select all entries.
% An omitted variable selects the first variable of the first selected basis.
% Explicit page and column order, including repeated indices, is preserved.
% Endpoint columns may be selected explicitly; they are never treated as an
% automatically accepted or truncated modal prefix.
%
% First derivatives use the existing solved-variable derivative evaluator.
% Boundary F derivatives use $$F_z=-N^2G/g$$. Other diagnostic derivatives
% require a separately implemented mathematical relation and are rejected.
%
% - Topic: Evaluate collections
% - Parameter z: physical coordinate column
% - Parameter options.variable: scalar variable name
% - Parameter options.derivativeOrder: requested vertical derivative order
% - Parameter options.pages: requested-page indices, empty means all
% - Parameter options.columns: explicit scientific column indices
% - Returns values: evaluated continuous basis columns by requested page
arguments (Input)
    self (1,1) IMBasisCollection
    z (:,1) double {mustBeReal,mustBeFinite}
    options.variable (1,1) string = ""
    options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.pages (1,:) double {mustBeInteger,mustBePositive} = []
    options.columns (1,:) double {mustBeInteger,mustBePositive} = []
end
arguments (Output)
    values double
end
pages = options.pages;
if isempty(pages)
    pages = 1:numel(self.basisIndex);
elseif any(pages > numel(self.basisIndex))
    error("IMBasisCollection:InvalidPage","pages must select existing requested pages.");
end
variable = options.variable;
if variable == ""
    variable = self.metadata(self.basisIndex(pages(1))).variables(1);
end
columnCounts = arrayfun(@(i) numel(self.metadata(self.basisIndex(i)).columnLabels),pages);
if isempty(options.columns)
    if any(columnCounts ~= columnCounts(1))
        error("IMBasisCollection:HeterogeneousColumns","Selected pages have different column counts. Evaluate homogeneous page groups or select explicit columns.");
    end
    columns = 1:columnCounts(1);
else
    columns = options.columns;
    if any(columns > min(columnCounts))
        error("IMBasisCollection:InvalidColumn","columns must exist on every selected page. Evaluate differing column counts separately.");
    end
end
values = zeros(numel(z),numel(columns),numel(pages));
% Evaluate each stored basis once per call, even for repeated requested pages.
for iBasis = unique(self.basisIndex(pages),"stable")
    basis = self.bases{iBasis};
    entry = self.metadata(iBasis);
    iVariable = find(entry.variables == variable,1);
    if isempty(iVariable)
        error("IMBasisCollection:UnsupportedVariable","Basis %d does not support %s. Inspect metadata.variables or evaluate compatible pages separately.",iBasis,variable);
    end
    if ~ismember(options.derivativeOrder,entry.derivativeOrders{iVariable})
        error("IMBasisCollection:UnsupportedDerivative","Basis %d does not implement derivative order %d for %s. Inspect metadata.derivativeOrders.",iBasis,options.derivativeOrder,variable);
    end
    if any(z < entry.zDomain(1) | z > entry.zDomain(2))
        error("IMBasisCollection:OutsideDomain","Evaluation coordinates must lie inside every selected basis domain.");
    end
    if options.derivativeOrder == 0
        sampled = basis.(variable)(z);
    elseif entry.columnKind == "endpoint"
        sampled = -(basis.N2(z)./basis.g).*basis.G(z);
    else
        sampled = basis.uz(z);
    end
    if size(sampled,1) ~= numel(z) || size(sampled,2) ~= numel(entry.columnLabels)
        error("IMBasisCollection:InvalidEvaluationShape","Basis %d evaluation must return nZ rows and one column per scientific label.",iBasis);
    end
    selected = find(self.basisIndex(pages) == iBasis);
    for iSelected = selected
        values(:,:,iSelected) = sampled(:,columns,self.sourcePage(pages(iSelected)));
    end
end
end
