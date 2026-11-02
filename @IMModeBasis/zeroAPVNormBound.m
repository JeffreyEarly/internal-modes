function bound = zeroAPVNormBound(self,boundary)
% Bound cancellation in canonical source form arithmetic, without a floor.
%
% The physical majorant may itself vanish in a neutral-label direction.
% Entrywise absolute canonical forms expose cancellation that a rounded
% congruence can otherwise report as a tiny positive norm.
rotation = self.zeroAPVSourceRotation_;
source = self.source.energyMatrix;
canonical = rotation'\source(:,:,1)/rotation;
absoluteForm = abs(canonical);
if isfinite(boundary.g0)
    source = self.source.surfaceBuoyancyMatrix;
    canonical = rotation'\source(:,:,1)/rotation;
    absoluteForm = absoluteForm+abs(boundary.g0)*abs(canonical);
end
if isfinite(boundary.gd)
    source = self.source.bottomBuoyancyMatrix;
    canonical = rotation'\source(:,:,1)/rotation;
    absoluteForm = absoluteForm+abs(boundary.gd)*abs(canonical);
end
absoluteMap = abs(self.zeroAPVColumnMap_);
bound = sum(absoluteMap.*(absoluteForm*absoluteMap),1);
end
