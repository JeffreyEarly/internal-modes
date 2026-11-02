function matrix = zeroAPVForm(self,name,isRaw)
% Pull an existing source form into current or unscaled profile coordinates.
if nargin < 3, isRaw = false; end
if self.family ~= "geostrophicZeroAPV"
    error("IMModeBasis:UnavailableQuadraticForm","Zero-APV quadratic forms require the geostrophicZeroAPV family.");
end
if self.nColumns == 0
    matrix = zeros(0);
    return
end
map = self.zeroAPVColumnMap_;
if ~isRaw, map = map./self.scaleFactors_; end
relativeMap = self.zeroAPVSourceRotation_\map;
sourceForm = self.source.(name);
matrix = relativeMap'*sourceForm(:,:,1)*relativeMap;
if any(~isfinite(matrix),"all")
    error("IMModeBasis:UnrepresentableQuadraticForm","The requested zero-APV form is not finite in these coordinates. Choose a representable normalization or pencil.");
end
matrix = IMGeostrophicZeroAPVFormTools.symmetrize(matrix);
end
