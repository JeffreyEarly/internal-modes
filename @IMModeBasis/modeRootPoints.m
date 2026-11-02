function [z,design] = modeRootPoints(self,generatingMode)
% Obtain physical endpoints and G roots from an already solved hydrostatic family.
%
% This generates points only. A separate weight fit completes an
% `IMVerticalGrid`. The requested internal mode must already be present;
% no auxiliary EVP is solved and no retained-band decision is made here.
%
% - Topic: Project modes
% - Parameter generatingMode: positive physical internal-mode number
% - Returns z: increasing physical sample points including both endpoints
% - Returns design: root-construction provenance
arguments
    self (1,1) IMModeBasis
    generatingMode (1,1) double {mustBeInteger,mustBePositive}
end
if self.family ~= "hydrostaticF" || ~self.includesBarotropic || ~isa(self.source,"IMBasisSet") || ~isa(self.source.solver,"IMSolverSpectral")
    error("IMModeBasis:UnavailableModeRoots","Mode-root points currently require a spectral rigid-lid hydrostatic F/G family.");
end
column = find(self.modeInfo.branch == "positive" & self.modeInfo.branchIndex == generatingMode,1);
if isempty(column)
    error("IMModeBasis:MissingGeneratingMode","Retain internal mode %d in the continuous family before requesting its roots.",generatingMode);
end
sourceColumn = self.columns(column)-1;
roots = sort(self.source.solver.rootsOfNativeMode(self.source.nativeModes(:,sourceColumn)));
tolerance = max(100*eps(max(1,max(abs(self.zDomain)))),1e-10*max(1,diff(self.zDomain)));
roots = roots(roots > self.zDomain(1)+tolerance & roots < self.zDomain(2)-tolerance);
roots = roots(:);
if ~isempty(roots), roots = roots([true;diff(roots) > tolerance]); end
z = [self.zDomain(1);roots;self.zDomain(2)];
design = IMDiscreteTransformTools.modeRootGridDesign(self.source,z,sourceColumn-1,generatingMode);
end
