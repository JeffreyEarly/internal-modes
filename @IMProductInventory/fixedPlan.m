function plan = fixedPlan(self,options)
% Reserve a deterministic sparse plan before numerical preparation.
%
% Supplied interaction IDs already obey the caller's physical rules. For
% each prefix, retained input factors use low/middle/cutoff ordinal stresses;
% their cumulative union preserves earlier checks. Fixed factors retain every
% column, including every boundary endpoint and frequency-sign coordinate.
%
% - Topic: Plan assessments
% - Parameter options.interactionIds: caller-selected valid IDs; default all
% - Parameter options.prefixCounts: increasing requested counts
% - Parameter options.productBudget: explicit reservation including zeros
% - Returns plan: metadata-only executable reservation
arguments (Input)
    self (1,1) IMProductInventory
    options.interactionIds (1,:) string = unique(self.products.interactionId,"stable").'
    options.prefixCounts (1,:) double {mustBeInteger,mustBePositive}
    options.productBudget (1,1) double {mustBeInteger,mustBePositive,mustBeFinite}
end
plan = IMProductPlan(self,options.interactionIds,options.prefixCounts,options.productBudget,"fixed");
end
