function transform = discreteTransform(self,grid)
% Build aligned sampled Galerkin operators from a continuous hydrostatic family.
%
% This fixed-grid bridge evaluates the existing basis without another EVP
% solve or weight fit. It preserves its labels and normalization. The
% coefficient map is weighted least squares; physical source pairings remain
% available separately through `modeProjectionFunctional`.
%
% The current bridge supports the rigid-lid hydrostatic F/G family. It is
% not a generic projection implementation for the other continuous families.
%
% - Topic: Project modes
% - Parameter grid: IMVerticalGrid with final physical points and weights
% - Returns transform: aligned IMInternalModesDiscreteTransform owning this basis
arguments
    self (1,1) IMModeBasis
    grid (1,1) IMVerticalGrid
end
if self.family ~= "hydrostaticF" || self.background.surfaceCondition ~= "rigidLid"
    error("IMModeBasis:UnavailableDiscreteTransform","The common fixed-grid transform bridge currently supports rigid-lid hydrostatic F/G modes only.");
end
if ~isequal(grid.zDomain,self.zDomain)
    error("IMModeBasis:IncompatibleDomain","The basis and sampling grid must have the same zDomain.");
end
if any(~ismember(self.modeInfo.branch,["barotropic","positive"]))
    error("IMModeBasis:InvalidHydrostaticFamily","A sampled hydrostatic family must select barotropic or internal columns.");
end
F = self.F(grid.z); G = self.G(grid.z);
[HF,HG] = self.hydrostaticNorms();
activeG = isfinite(self.h);
N2 = self.background.N2(grid.z);
interiorF = ones(size(grid.z)); interiorG = N2/self.background.g;
channel = @(interior,target,active) struct(available=true,reason="",activeModeMask=active, ...
    interiorWeight=interior,endpointMetricMatrix=zeros(grid.nPoints),metricMatrix=diag(interior.*grid.weights), ...
    targetGramMatrix=diag(target),targetMajorantGramMatrix=diag(target));
channels = struct(F=channel(interiorF,HF,true(1,self.nColumns)),G=channel(interiorG,HG,activeG));
metadata = struct(construction="IMModeBasis.discreteTransform",background=self.background,surfaceCondition="rigidLid", ...
    columnLabels=self.columnLabels,relation="F = h dG/dz");
transform = IMInternalModesDiscreteTransform(z=grid.z,weights=grid.weights,modeNumber=self.modeInfo.branchIndex.',h=self.h, ...
    normalization=self.normalization,inverseF=F,inverseG=G,endpointF=self.F(self.zDomain([2 1]).'),endpointG=self.G(self.zDomain([2 1]).'), ...
    channelData=channels,primaryVariable="G",zDomain=self.zDomain,g=self.background.g,modeFamily="hydrostatic", ...
    N2Values=N2,problemMetadata=metadata,basis=self);
end
