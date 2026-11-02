function basis = rotateSurfaceBuoyancy(self,options)
% Diagonalize surface buoyancy against signed generalized energy.
%
% The current retained pencil is g B0 c = chi Hg c, ordered with nonzero
% surface response first. The right norms have unit absolute value and
% preserve their signs. The canonical source must include a surface
% coordinate, including when the retained inventory is empty.
%
% - Topic: Rotate zero-APV modes
% - Parameter options.g0: optional signed surface coefficient; retained specification by default
% - Parameter options.gd: optional signed bottom coefficient; retained specification by default
% - Returns basis: immutable surface-buoyancy coordinates on the same span
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
if ~ismember("surface",self.canonicalEndpointLabels)
    error("IMModeBasis:SurfaceEndpointRequired","Surface-buoyancy rotation requires a surface endpoint source coordinate.");
end
right = self.generalizedEnergyMatrix(g0=boundary.g0,gd=boundary.gd);
basis = self.rotateZeroAPVPencil("surfaceBuoyancy",self.background.g*self.surfaceBuoyancyMatrix,right,"nonzeroFirst","unitSignedGeneralizedEnergy",boundary.g0,boundary.gd);
end
