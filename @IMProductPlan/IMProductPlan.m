classdef IMProductPlan
    % Reserve individual product identities before invoking expensive callbacks.
    %
    % The stored inventory is immutable. Execution caches are local to assess;
    % a later execution cannot reuse stale evaluated fields or reference data.
    %
    % ```matlab
    % plan = inventory.fixedPlan(prefixCounts=1:8,productBudget=10000);
    % result = plan.assess(grids,chunkSize=128);
    % ```
    %
    % - Topic: Inspect reservations
    % - Topic: Execute assessments
    % - Declaration: classdef IMProductPlan
    properties (SetAccess = private)
        % Original metadata-only inventory.
        % - Topic: Inspect reservations
        inventory
        % Selection rule, fixed sparse or bounded all-products validation.
        % - Topic: Inspect reservations
        selection
        % Examined retained-count values, unchanged by assessment decisions.
        % - Topic: Inspect reservations
        prefixCounts
        % Product-family rows selected in original inventory order.
        % - Topic: Inspect reservations
        rows
        % Two-by-product input column indices for each selected row.
        % - Topic: Inspect reservations
        pairs
        % Earliest retained count at which each selected product is examined.
        % - Topic: Inspect reservations
        firstCounts
        % Exact reserved product count, including future structural zeros.
        % - Topic: Inspect reservations
        reservedProducts
        % Explicit maximum authorized product count.
        % - Topic: Inspect reservations
        productBudget
        % Inventory construction/planning elapsed seconds.
        % - Topic: Inspect reservations
        planningSeconds
    end
    methods
        function self = IMProductPlan(inventory,interactionIds,prefixCounts,productBudget,selection)
            % Validate and reserve a complete metadata-only execution plan.
            % - Topic: Inspect reservations
            arguments (Input)
                inventory (1,1) IMProductInventory
                interactionIds (1,:) string
                prefixCounts (1,:) double {mustBeInteger,mustBePositive}
                productBudget (1,1) double {mustBeInteger,mustBePositive,mustBeFinite}
                selection (1,1) string {mustBeMember(selection,["fixed","all"])}
            end
            timer = tic;
            if isempty(prefixCounts) || any(diff(prefixCounts)<=0) || isempty(interactionIds) || numel(unique(interactionIds))~=numel(interactionIds) || any(~ismember(interactionIds,inventory.products.interactionId))
                error("IMProductPlan:InvalidSelection","Counts must strictly increase; interaction IDs must be distinct supplied identities.");
            end
            rows = find(ismember(inventory.products.interactionId,interactionIds)).';
            pairs = cell(size(rows)); firstCounts = cell(size(rows)); reserved = 0;
            for j=1:numel(rows)
                row = inventory.products(rows(j),:);
                a = inventory.factors{row.factorA}; b = inventory.factors{row.factorB};
                for factor={a,b}
                    if factor{1}.countRole=="retained" && ~all(ismember(1:prefixCounts(end),factor{1}.ordinals))
                        error("IMProductPlan:InsufficientInputBand","Retained factors must contain every explicitly requested input ordinal; no count is silently reduced.");
                    end
                end
                [pairs{j},firstCounts{j},number] = reservePairs(a,b,prefixCounts(end),selection,productBudget-reserved);
                reserved=reserved+number;
                output = inventory.outputs{row.output};
                if output.countRole=="retained" && ~all(ismember(1:prefixCounts(end),output.ordinals))
                    error("IMProductPlan:InvalidOutputCounts","Every ordinal through the requested count must exist in each retained output; sparse count lists do not permit missing output modes.");
                end
            end
            if reserved==0
                error("IMProductPlan:EmptySelection","The selected counts contain no input products.");
            end
            self.inventory=inventory; self.selection=selection; self.prefixCounts=prefixCounts;
            self.rows=rows; self.pairs=pairs; self.firstCounts=firstCounts;
            self.reservedProducts=reserved; self.productBudget=productBudget; self.planningSeconds=toc(timer);
        end
    end
end

function [pairs,firstCounts,number] = reservePairs(a,b,maximum,selection,budget)
% Enumerate small ordinal stresses first; expand sign/column multiplicities
% only after their reservation fits. Large unrequested basis tails cost no
% Cartesian product storage.
if selection=="all"
    ia=1:numel(a.labels); ib=1:numel(b.labels);
    if a.countRole=="retained", ia=ia(a.ordinals<=maximum); end
    if b.countRole=="retained", ib=ib(b.ordinals<=maximum); end
    number=numel(ia)*numel(ib); requireBudget(number,budget);
    [i,j]=ndgrid(ia,ib); pairs=[i(:).';j(:).']; firstCounts=ones(1,number);
    if a.countRole=="retained", firstCounts=max(firstCounts,a.ordinals(pairs(1,:))); end
    if b.countRole=="retained", firstCounts=max(firstCounts,b.ordinals(pairs(2,:))); end
    return
end
ordinalPairs=zeros(0,2); first=zeros(0,1); multiplicities=zeros(0,1); number=0;
for count=1:maximum
    anchors=unique([1 2 ceil(count/2) max(1,count-2):count]); anchors=anchors(anchors<=count);
    left=anchors; right=anchors;
    if a.countRole=="fixed", left=0; end
    if b.countRole=="fixed", right=0; end
    [i,j]=ndgrid(left,right); candidates=[i(:),j(:)];
    additions=candidates(~ismember(candidates,ordinalPairs,"rows"),:);
    for q=1:size(additions,1)
        na=nnz(a.ordinals==additions(q,1)); nb=nnz(b.ordinals==additions(q,2));
        if additions(q,1)==0, na=numel(a.labels); end
        if additions(q,2)==0, nb=numel(b.labels); end
        amount=na*nb; number=number+amount; requireBudget(number,budget);
        ordinalPairs(end+1,:)=additions(q,:); first(end+1,1)=count; multiplicities(end+1,1)=amount; %#ok<AGROW>
    end
end
pairs=zeros(2,number); firstCounts=zeros(1,number); position=0;
for q=1:size(ordinalPairs,1)
    ia=find(a.ordinals==ordinalPairs(q,1)); ib=find(b.ordinals==ordinalPairs(q,2));
    if ordinalPairs(q,1)==0, ia=1:numel(a.labels); end
    if ordinalPairs(q,2)==0, ib=1:numel(b.labels); end
    [i,j]=ndgrid(ia,ib); batch=position+(1:multiplicities(q));
    pairs(:,batch)=[i(:).';j(:).']; firstCounts(batch)=first(q); position=position+multiplicities(q);
end
[~,order]=sortrows(pairs.',[2 1]); pairs=pairs(:,order); firstCounts=firstCounts(order);
end

function requireBudget(number,budget)
if number>budget
    error("IMProductPlan:ProductBudgetExceeded","The metadata reservation exceeds the remaining product budget of %d (at least %d products, including zeros). No evaluator has run.",budget,number);
end
end
