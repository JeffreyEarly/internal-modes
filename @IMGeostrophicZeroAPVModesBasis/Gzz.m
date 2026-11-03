function values = Gzz(self,z,options)
% Evaluate the regular diagnostic derivative Gzz.
% - Topic: Evaluate geostrophic zero-APV modes
arguments
    self IMGeostrophicZeroAPVModesBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    options.pages (1,:) double {mustBeInteger,mustBePositive} = 1:numel(self.k)
end
values = self.evaluate(z,variable="G",derivativeOrder=2,pages=options.pages);
end
