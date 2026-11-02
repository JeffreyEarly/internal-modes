function basis = rotateBoundaryDepth(self,options)
% Diagonalize signed generalized energy against endpoint response.
%
% The pencil uses only the current retained span. Its eigenvalues gamma
% define h0 = 2 k^2 gamma, distinct from ordinary EVP equivalent depths.
% The output has unit endpoint-response right norms before any later
% positive profile normalization.
%
% - Topic: Rotate zero-APV modes
% - Parameter options.g0: optional signed surface coefficient; retained specification by default
% - Parameter options.gd: optional signed bottom coefficient; retained specification by default
% - Returns basis: immutable boundary-depth coordinates on the same span
arguments (Input)
    self (1,1) IMModeBasis
    options.g0 double {mustBeReal} = []
    options.gd double {mustBeReal} = []
end
arguments (Output)
    basis (1,1) IMModeBasis
end
if self.family ~= "geostrophicZeroAPV"
    error("IMModeBasis:UnavailableRotation","Zero-APV rotations require the geostrophicZeroAPV family.");
end
boundary = self.zeroAPVEnergyBoundary(options.g0,options.gd);
left = self.generalizedEnergyMatrix(g0=boundary.g0,gd=boundary.gd);
basis = self.rotateZeroAPVPencil("boundaryDepth",left,self.endpointResponseMetric,"descending","unitEndpointResponseMetric",boundary.g0,boundary.gd);
end
