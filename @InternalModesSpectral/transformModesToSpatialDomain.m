function [F,G] = transformModesToSpatialDomain(self,V_cheb,h,maxModes)

F = zeros(length(self.z),maxModes);
G = zeros(length(self.z),maxModes);
h = reshape(h(1:maxModes),1,[]);

varargout = cell(size(varargin));
for iArg=1:length(varargin)
    varargout{iArg} = zeros(1,maxModes);
end

% This still need to be optimized to *not* do the transforms
% twice, when the EVP grid is the same as the output grid.
[maxIndexZ] = find(self.N2_xLobatto-self.gridFrequency*self.gridFrequency>0,1,'first');
if maxIndexZ > 1 % grab a point just above the turning point, which should have the right sign.
    maxIndexZ = maxIndexZ-1;
elseif isempty(maxIndexZ)
    maxIndexZ = 1;
end
for j=1:maxModes
    Fj = self.FFromVCheb(V_cheb(:,j),h(j));
    Gj = self.GFromVCheb(V_cheb(:,j),h(j));
    switch self.normalization
        case Normalization.uMax
            A = max( abs( Fj ));
        case Normalization.wMax
            A = max( abs( Gj ) );
        case Normalization.kConstant
            A = sqrt(self.GNorm( Gj ));
        case Normalization.omegaConstant
            A = sqrt(self.FNorm( Fj ));
        case Normalization.geostrophic
            A = sqrt(self.GeostrophicNorm( Gj ));
    end
    if Fj(maxIndexZ) < 0
        A = -A;
    end

    G(:,j) = self.GOutFromVCheb(V_cheb(:,j),h(j))/A;
    F(:,j) = self.FOutFromVCheb(V_cheb(:,j),h(j))/A;
end

end