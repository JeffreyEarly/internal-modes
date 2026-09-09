function plan = allProductsPlan(self,options)
% Reserve all supplied products for bounded independent validation.
%
% This validation control measures the supplied inventory only. It adds no
% physical interactions, universal coverage, or superposition guarantee.
%
% - Topic: Plan assessments
% - Parameter options.interactionIds: independently chosen supplied interactions
% - Parameter options.prefixCounts: increasing requested counts
% - Parameter options.productBudget: explicit reservation including zeros
% - Returns plan: bounded validation reservation
arguments (Input)
    self (1,1) IMProductInventory
    options.interactionIds (1,:) string = unique(self.products.interactionId,"stable").'
    options.prefixCounts (1,:) double {mustBeInteger,mustBePositive}
    options.productBudget (1,1) double {mustBeInteger,mustBePositive,mustBeFinite}
end
plan = IMProductPlan(self,options.interactionIds,options.prefixCounts,options.productBudget,"all");
end
