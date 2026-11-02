function matrix = generalizedEnergyMatrix(self,options)
% Form signed generalized energy in the current retained coordinates.
%
% Defaults retain the scientific g0/gd coefficients. Absent canonical
% endpoints contribute no term, including when their coefficient is Inf.
% This operation cannot remove or create a source endpoint coordinate.
%
% - Topic: Form zero-APV quadratic forms
% - Parameter options.g0: optional signed surface coefficient; defaults to the retained specification
% - Parameter options.gd: optional signed bottom coefficient; defaults to the retained specification
% - Returns matrix: H + g0 B0 + gd Bd on the current retained span
arguments (Input)
    self (1,1) IMModeBasis
    options.g0 double {mustBeReal} = []
    options.gd double {mustBeReal} = []
end
arguments (Output)
    matrix (:,:) double
end
boundary = self.zeroAPVEnergyBoundary(options.g0,options.gd);
matrix = self.energyMatrix;
if isfinite(boundary.g0) && ismember("surface",self.canonicalEndpointLabels)
    matrix = matrix+boundary.g0*self.surfaceBuoyancyMatrix;
end
if isfinite(boundary.gd) && ismember("bottom",self.canonicalEndpointLabels)
    matrix = matrix+boundary.gd*self.bottomBuoyancyMatrix;
end
if any(~isfinite(matrix),"all")
    error("IMModeBasis:UnrepresentableQuadraticForm","Generalized energy is not finite in these coefficient coordinates.");
end
matrix = IMGeostrophicZeroAPVFormTools.symmetrize(matrix);
end
