function info = samplingIntegrationInfo(self)
% Describe immutable source resolution and known integration subdivisions.
%
% The original inventory and actual retained source survive column selection,
% normalization and zero-APV rotation. They seed quadrature effort without
% asserting convergence of arbitrary profiles or intrinsic evaluators.
% - Topic: Developer topics
% - Returns info: source inventory/resolution and known physical breakpoints
% - Developer: true
arguments (Input)
    self (1,1) IMModeBasis
end
arguments (Output)
    info (1,1) struct
end
points = self.integrationIntervals_(:);
resolved = self.scientificSpecification.resolved;
if isfield(resolved,"intervals"), points = [points;resolved.intervals(:)]; end
sourceClass = "intrinsic";
sourceResolution = 0;
if ~isempty(self.source)
    sourceClass = string(class(self.source));
    if isprop(self.source,"solver")
        if isprop(self.source.solver,"nEVP"), sourceResolution = self.source.solver.nEVP; end
        if isa(self.source.solver,"IMSolverFiniteDifference"), points = [points;self.source.solver.zNative]; end
    end
    if isa(self.source,"IMGeostrophicModesBasis") || isa(self.source,"IMGeostrophicZeroAPVModesBasis")
        points = [points;self.source.samplingBreakpoints()];
    end
end
points = unique(points(points > self.zDomain(1) & points < self.zDomain(2)));
info = struct(zDomain=self.zDomain,breakpoints=points,originalColumnCount=self.normalizationColumnCount_,sourceResolution=sourceResolution,sourceClass=sourceClass,resolutionStatus="integration effort hint; not a convergence qualification");
end
