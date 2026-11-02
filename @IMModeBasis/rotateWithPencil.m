function basis = rotateWithPencil(self,options)
% Rotate the current retained span with a real symmetric matrix pencil.
%
% Both matrices use current normalized coefficient coordinates. The
% operation never restores discarded modes. The output retains signed
% right-unit normalization, which can be restored with rotationRightMetric
% after an explicit positive profile normalization.
%
% - Topic: Rotate zero-APV modes
% - Parameter options.name: nonempty scientific name for this rotation
% - Parameter options.leftMatrix: current-column symmetric left form
% - Parameter options.rightMatrix: current-column symmetric right form
% - Returns basis: immutable rotated common basis
arguments (Input)
    self (1,1) IMModeBasis
    options.name (1,1) string
    options.leftMatrix (:,:) double {mustBeReal,mustBeFinite}
    options.rightMatrix (:,:) double {mustBeReal,mustBeFinite}
end
arguments (Output)
    basis (1,1) IMModeBasis
end
if ismissing(options.name) || strlength(strtrim(options.name)) == 0
    error("IMModeBasis:InvalidRotationName","Supply a nonempty rotation name.");
end
ordering = "descending";
if self.family == "geostrophicZeroAPV" && IMGeostrophicZeroAPVFormTools.pagesEqual(options.leftMatrix,self.background.g*self.surfaceBuoyancyMatrix)
    ordering = "nonzeroFirst";
end
basis = self.rotateZeroAPVPencil(options.name,options.leftMatrix,options.rightMatrix,ordering,"unitSignedRightMetric",[],[]);
end
