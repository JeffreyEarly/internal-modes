function [transform,assessment] = discreteTransform(self,varargin)
% Build a fixed-grid transform or use the inherited qualification API.
%
% A positional IMVerticalGrid uses its points and weights exactly,
% without a solve or weight fit. Name-value requests retain the
% parent's full transform qualification interface.
%
% - Topic: Build discrete transforms
% - Parameter varargin: one IMVerticalGrid or inherited name-value options
% - Returns transform: aligned F/G transform with inactive G0
% - Returns assessment: inherited qualification result, empty for a fixed grid
if isempty(varargin) || ~isa(varargin{1},"IMVerticalGrid")
    [transform,assessment] = discreteTransform@IMInternalModesBasis(self,varargin{:});
    return
end
if ~isscalar(varargin{1})
    error("IMHydrostaticModesBasis:InvalidGrid","Supply one scalar IMVerticalGrid.");
end
grid = varargin{1};
options = parseRecoveryOptions(varargin{2:end});
if ~isequal(grid.zDomain,self.zDomain)
    error("IMHydrostaticModesBasis:IncompatibleDomain","The basis and sampling grid must have the same zDomain.");
end
F = self.F(grid.z); G = self.G(grid.z);
if ismember(string(self.normalization),["depth","depthRMSF"])
    HF = diff(self.zDomain)*ones(1,numel(self.h));
else
    HF = diag(self.gramMatrix(variable="F")).';
end
HG = HF./self.h;
HG(~isfinite(self.h)) = 0;
activeG = isfinite(self.h);
N2 = self.N2(grid.z);
channel = @(interior,target,active) struct(available=true,reason="",activeModeMask=active,interiorWeight=interior,endpointMetricMatrix=zeros(grid.nPoints),metricMatrix=diag(interior.*grid.weights),targetGramMatrix=diag(target),targetMajorantGramMatrix=diag(target));
channels = struct(F=channel(ones(size(grid.z)),HF,true(1,numel(self.h))),G=channel(N2/self.evp.g,HG,activeG));
if ~isempty(options.recoveryF), channels.F.recovery = options.recoveryF; end
if ~isempty(options.recoveryG), channels.G.recovery = options.recoveryG; end
if options.forwardOperator ~= ""
channels.F.forwardOperator = options.forwardOperator;
channels.G.forwardOperator = options.forwardOperator;
end
metadata = struct(construction="IMHydrostaticModesBasis.discreteTransform",surfaceCondition="rigidLid",columnLabels=self.columnLabels,relation="F = h dG/dz with F0=1 and G0=0");
transform = IMInternalModesDiscreteTransform(z=grid.z,weights=grid.weights,modeNumber=self.modeNumber,h=self.h,normalization=string(self.normalization),inverseF=F,inverseG=G,endpointF=self.F(self.zDomain([2 1]).'),endpointG=self.G(self.zDomain([2 1]).'),channelData=channels,primaryVariable="G",zDomain=self.zDomain,g=self.evp.g,modeFamily="hydrostatic",N2Values=N2,problemMetadata=metadata,basis=self);
assessment = [];
end

function options = parseRecoveryOptions(options)
arguments (Input)
    options.recoveryF = []
    options.recoveryG = []
    options.forwardOperator (1,1) string = ""
end
end
