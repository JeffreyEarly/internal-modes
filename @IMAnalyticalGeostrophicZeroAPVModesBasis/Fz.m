function values = Fz(self,z,options)
% Evaluate the regular diagnostic derivative Fz.
% - Topic: Evaluate geostrophic zero-APV modes
arguments
    self IMAnalyticalGeostrophicZeroAPVModesBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    options.pages (1,:) double {mustBeInteger,mustBePositive} = 1:numel(self.k)
end
values = self.evaluate(z,variable="F",derivativeOrder=1,pages=options.pages);
end
