function [projection,recipe] = prepareProjection(self,z,weights,options)
% Construct a fixed projection for one requested page and variable.
%
% This operation preserves the supplied samples, weights, and explicit
% columns. It performs no quadrature fitting or acceptance policy. The
% continuous basis supplies the signed pairing and positive error metric.
%
% - Topic: Construct projections
% - Declaration: [projection,recipe] = prepareProjection(collection,z,weights,options)
% - Parameter z: physical sample column
% - Parameter weights: fixed quadrature weights aligned with z
% - Parameter options.page: requested-page index, not distinct-solve index
% - Parameter options.variable: scalar u or aligned F/G variable
% - Parameter options.columns: explicit column indices; all by default, [] selects none
% - Parameter options.prepareRecovery: also prepare numerical recovery for Galerkin records
% - Returns projection: IMProjection with the requested columns
% - Returns recipe: basis-owned metric, capabilities, and provenance
arguments (Input)
    self (1,1) IMBasisCollection
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite}
    options.page (1,1) double {mustBeInteger,mustBePositive} = 1
    options.variable (1,1) string = ""
    options.columns = ":"
    options.prepareRecovery (1,1) logical = true
end
if options.page > numel(self.basisIndex)
    error("IMBasisCollection:InvalidPage","page must select an existing requested page.");
end
columns = options.columns;
if isequal(columns,":") || isequal(columns,':')
    columns = self.allColumnIndices(options.page);
end
if ~isa(columns,"double") || ~isrow(columns) && ~isempty(columns) || ~isreal(columns) ...
        || any(~isfinite(columns) | columns < 1 | columns ~= fix(columns),"all")
    error("IMProjection:InvalidColumns","columns must be a row of positive integer indices, [] or ':' for all columns.");
end
columns = reshape(columns,1,[]);
basis = self.bases{self.basisIndex(options.page)}.selectColumns(columns);
if options.variable == ""
    recipe = basis.projectionRecipe();
else
    recipe = basis.projectionRecipe(variable=options.variable);
end
if isa(basis,"IMNeutralLabelModesBasis")
    projection = basis.projectionOnGrid(z,weights,variable="APV",prepareRecovery=options.prepareRecovery);
else
    % Auxiliary EVP recipes remain useful for Gram qualification. They do
    % not imply a complete physical state/source projector for that family.
    projection = IMProjection.fromRecipe(recipe,z,weights,prepareRecovery=options.prepareRecovery);
end
end
