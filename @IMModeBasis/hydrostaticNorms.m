function [HF,HG] = hydrostaticNorms(self)
% Compute diagonal continuous F/G norms in the captured common normalization.
arguments
    self (1,1) IMModeBasis
end
if self.family ~= "hydrostaticF" || self.background.surfaceCondition ~= "rigidLid"
    error("IMModeBasis:UnavailableDiscreteTransform","The common F/G metric currently supports rigid-lid hydrostatic modes only.");
end
if self.normalization == "depthRMSF"
    HF = self.background.D*ones(1,self.nColumns);
else
    [z,w] = self.normalizationQuadrature(false);
    HF = sum(w.*self.F(z).^2,1);
end
HG = HF./self.h;
HG(~isfinite(self.h)) = 0;
end
