function points = samplingBreakpoints(self)
% Return the actual interpolation knots of both finite-difference fields.
%
% These subdivisions describe retained source interpolation, independently
% of the selected or rotated endpoint columns and observation geometry.
% - Topic: Developer topics
% - Returns points: sorted native F/G finite-difference knots
% - Developer: true
arguments (Input)
    self (1,1) IMGeostrophicZeroAPVModesBasis
end
arguments (Output)
    points (:,1) double
end
points = zeros(0,1);
if isa(self.solver,"IMSolverFiniteDifference"), points = self.solver.zNative; end
if isa(self.gSolver,"IMSolverFiniteDifference"), points = [points;self.gSolver.zNative]; end
points = unique(points);
end
