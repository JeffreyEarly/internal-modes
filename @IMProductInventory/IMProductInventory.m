classdef IMProductInventory
    % Declare products, continuous factors, and caller-owned output recipes.
    %
    % Construction inspects metadata only. Physical channels, valid
    % interactions, frequency signs, and endpoint recipes belong to the
    % caller. Factors and output recipes are evaluated only after a plan
    % reserves every selected individual product, including exact zeros.
    %
    % ```matlab
    % inventory = IMProductInventory(factors,products,outputs);
    % plan = inventory.fixedPlan(prefixCounts=1:8,productBudget=10000);
    % ```
    %
    % - Topic: Declare inventories
    % - Topic: Plan assessments
    % - Declaration: classdef IMProductInventory
    properties (SetAccess = private)
        % Cell row of factor metadata and explicit evaluation callbacks.
        % - Topic: Declare inventories
        factors
        % Ordered table of valid interaction/channel/factor/output identities.
        % - Topic: Declare inventories
        products
        % Cell row of metadata and deferred output preparation callbacks.
        % - Topic: Declare inventories
        outputs
    end
    methods
        function self = IMProductInventory(factors,products,outputs)
            % Validate an inventory without invoking any evaluator.
            %
            % Factors require id, family, labels, ordinals, frequencySigns,
            % countRole (retained/fixed), evaluate(z,columns,referenceId),
            % and provenance. Outputs require id, family, labels, ordinals,
            % countRole, prepare(grids), and provenance. Product table columns
            % are interactionId, channel, factorA, factorB, output; numeric
            % references index the corresponding cell rows. Optional complex
            % coefficient multiplies the entire product on every grid.
            %
            % - Topic: Declare inventories
            % - Parameter factors: cell row of deferred factors
            % - Parameter products: ordered valid product-family table
            % - Parameter outputs: cell row of deferred output projections
            arguments (Input)
                factors (1,:) cell
                products table
                outputs (1,:) cell
            end
            if isempty(factors) || isempty(outputs) || isempty(products)
                error("IMProductInventory:EmptyInventory","Supply factors, outputs, and at least one valid product-family row.");
            end
            required = ["interactionId","channel","factorA","factorB","output"];
            if ~all(ismember(required,string(products.Properties.VariableNames)))
                error("IMProductInventory:InvalidProducts","Product table requires interactionId, channel, factorA, factorB, and output.");
            end
            for kind = ["factor","output"]
                if kind == "factor", entries = factors; else, entries = outputs; end
                ids = strings(1,numel(entries));
                for i = 1:numel(entries)
                    entry = entries{i};
                    fields = ["id","family","labels","ordinals","countRole","provenance"];
                    if kind == "factor", fields = [fields,"frequencySigns","evaluate"]; else, fields = [fields,"prepare"]; end %#ok<AGROW>
                    if ~isstruct(entry) || ~isscalar(entry) || ~all(isfield(entry,fields))
                        error("IMProductInventory:InvalidEntry","Each %s needs %s.",kind,join(fields,", "));
                    end
                    entry.id = string(entry.id); entry.family = string(entry.family);
                    entry.labels = reshape(string(entry.labels),1,[]);
                    entry.ordinals = reshape(entry.ordinals,1,[]);
                    entry.countRole = string(entry.countRole);
                    if ~isscalar(entry.id) || ismissing(entry.id) || strlength(entry.id)==0 || ~isscalar(entry.family) || ismissing(entry.family) || strlength(entry.family)==0 || isempty(entry.labels) || any(ismissing(entry.labels) | strlength(entry.labels)==0) || numel(unique(entry.labels)) ~= numel(entry.labels)
                        error("IMProductInventory:InvalidIdentity","Entry identities and distinct column labels must be nonempty.");
                    end
                    if ~isnumeric(entry.ordinals) || ~isreal(entry.ordinals) || numel(entry.ordinals) ~= numel(entry.labels) || any(~isfinite(entry.ordinals) | entry.ordinals < 1 | fix(entry.ordinals) ~= entry.ordinals) || ~isscalar(entry.countRole) || ~ismember(entry.countRole,["retained","fixed"])
                        error("IMProductInventory:InvalidOrdinals","Column ordinals must be positive integers; countRole must be retained or fixed.");
                    end
                    if ~isstruct(entry.provenance) || ~isscalar(entry.provenance) || isempty(fieldnames(entry.provenance))
                        error("IMProductInventory:MissingProvenance","Every factor/output requires explicit provenance.");
                    end
                    if isfield(entry,"columnKind") && string(entry.columnKind)=="endpoint" && entry.countRole~="fixed"
                        error("IMProductInventory:EndpointPrefix","Boundary endpoint coordinates require fixed countRole.");
                    end
                    if kind == "factor"
                        entry.frequencySigns = reshape(entry.frequencySigns,1,[]);
                        if ~isa(entry.evaluate,"function_handle") || numel(entry.frequencySigns) ~= numel(entry.labels) || ~isreal(entry.frequencySigns) || any(~ismember(entry.frequencySigns,[-1 0 1]))
                            error("IMProductInventory:InvalidFactor","Factors require an evaluator and one frequency sign per column.");
                        end
                        factors{i} = entry;
                    else
                        if ~isa(entry.prepare,"function_handle") || any(diff(entry.ordinals)<0)
                            error("IMProductInventory:InvalidOutput","Outputs require a preparation callback and nondecreasing ordinals.");
                        end
                        outputs{i} = entry;
                    end
                    ids(i) = entry.id;
                end
                if numel(unique(ids))~=numel(ids)
                    error("IMProductInventory:DuplicateIdentity","%s IDs must be distinct.",kind);
                end
            end
            products.interactionId = string(products.interactionId);
            products.channel = string(products.channel);
            if size(products.interactionId,2)~=1 || size(products.channel,2)~=1 || any(ismissing(products.interactionId) | strlength(products.interactionId)==0 | ismissing(products.channel) | strlength(products.channel)==0)
                error("IMProductInventory:InvalidProducts","Interaction and channel identities must be nonempty scalar strings per row.");
            end
            for field = ["factorA","factorB","output"]
                values = products.(field);
                bound = numel(factors); if field=="output", bound=numel(outputs); end
                if ~isnumeric(values) || ~isreal(values) || size(values,2)~=1 || any(~isfinite(values) | values<1 | values>bound | fix(values)~=values)
                    error("IMProductInventory:InvalidProducts","%s must index declared entries.",field);
                end
            end
            if ~ismember("coefficient",string(products.Properties.VariableNames))
                products.coefficient = ones(height(products),1);
            elseif ~isnumeric(products.coefficient) || size(products.coefficient,2)~=1 || any(~isfinite(products.coefficient))
                error("IMProductInventory:InvalidProducts","Product coefficients must be finite scalars per row.");
            end
            self.factors = factors; self.products = products; self.outputs = outputs;
        end
    end
    methods (Static)
        factor = collectionFactor(collection,options)
    end
end
