% Take matrices A and B from the generalized eigenvalue problem
% (GEP) and returns F,G,h. The last seven arguments are all
% function handles that do as they say.
function [F,G,h,varargout] = ModesFromGEP(self,A,B,varargin,options)
arguments
    self InternalModesSpectral
    A (:,:) double
    B (:,:) double
end
arguments (Repeating)
    varargin
end
arguments
    options.negativeEigenvalues = 0
end
if ( any(any(isnan(A))) || any(any(isnan(B))) )
    error('EVP setup fail. Found at least one nan in matrices A and B.\n');
end
[V,D] = eig( A, B );

% The following might be better, as it captures the the
% barotopic mode.
% d = diag(D);
% [d, permutation] = sort(real(d),'ascend');
% V_cheb=V(:,permutation);
% if options.negativeEigenvalues > 0
%     negIndices = find(d<0,options.negativeEigenvalues,'last');
% else
%     negIndices = [];
% end
% if isempty(negIndices)
%     minIndex = find(d>=0,1,'first');
%     if isempty(minIndex)
%         fprintf('No usable modes found! Try with higher resolution.\n');
%         return;
%     end
% else
%     minIndex = min(negIndices);
% end
% V_cheb = V_cheb(:,minIndex:end);
% h = self.hFromLambda(d(minIndex:end));

[h, permutation] = sort(real(self.hFromLambda(diag(D))),'descend');
V_cheb=V(:,permutation);
if options.negativeEigenvalues > 0
    negIndices = find(h<0,options.negativeEigenvalues,'first');
    permutation = cat(1,negIndices,setdiff((1:length(h))',negIndices));
    h = h(permutation);
    V_cheb=V_cheb(:,permutation);
end

if self.nModes == 0
    maxModes = ceil(find(h>0,1,'last')/2); % Have to do ceil, not floor, or we lose the barotropic mode.
    if maxModes == 0
        fprintf('No usable modes found! Try with higher resolution.\n');
        return;
    end
else
    maxModes = self.nModes;
end

%             FzOutFromGCheb = @(G_cheb,h) h * self.T_xCheb_zOut(self.Diff1_xCheb(self.Diff1_xCheb(G_cheb)));
%             Fz = zeros(length(self.z),maxModes);

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
    %                 Fz(:,j) = FzOutFromGCheb(G_cheb(:,j),h(j))/A;
    % K-constant norm: G(0)^2 + \frac{1}{g} \int_{-D}^0 (N^2 -
    % f_0^2)
    for iArg=1:length(varargin)
        if ( strcmp(varargin{iArg}, 'F2') )
            varargout{iArg}(j) = self.Lz*self.FNorm( Fj/A );
        elseif ( strcmp(varargin{iArg}, 'G2') )
            varargout{iArg}(j) = self.Lz*self.FNorm( Gj/A );
        elseif ( strcmp(varargin{iArg}, 'N2G2') )
            varargout{iArg}(j) = self.g*(self.GNorm( Gj/A )-Gj(1)*Gj(1)) + self.f0*self.f0*self.Lz*self.FNorm( Gj/A ); % this is being clever, but should give \int N2*G2 dz
        elseif  ( strcmp(varargin{iArg}, 'uMax') )
            B = max( abs( Fj ));
            varargout{iArg}(j) = abs(A/B);
        elseif  ( strcmp(varargin{iArg}, 'wMax') )
            B = max( abs( Gj ) );
            varargout{iArg}(j) = abs(A/B);
        elseif ( strcmp(varargin{iArg}, 'kConstant') )
            B = sqrt(self.GNorm( Gj ));
            varargout{iArg}(j) = abs(A/B);
        elseif ( strcmp(varargin{iArg}, 'omegaConstant') )
            B = sqrt(self.FNorm( Fj ));
            varargout{iArg}(j) = abs(A/B);
        elseif ( strcmp(varargin{iArg}, 'geostrophicNorm') )
            B = sqrt(self.GeostrophicNorm( Gj ));
            varargout{iArg}(j) = abs(A/B);
        elseif ( strcmp(varargin{iArg}, 'int_N2_G_dz/g') )
            varargout{iArg}(j) = sum(self.Int_xCheb .*InternalModesSpectral.fct((1/self.g) * self.N2_xLobatto .* (Gj/A)));
        elseif ( strcmp(varargin{iArg}, 'int_F_dz') )
            varargout{iArg}(j) = sum(self.Int_xCheb .*InternalModesSpectral.fct(Fj/A));
        elseif ( strcmp(varargin{iArg}, 'int_G_dz') )
            varargout{iArg}(j) = sum(self.Int_xCheb .*InternalModesSpectral.fct(Gj/A));
        else
            error('Invalid option. You may request F2, G2, N2G2');
        end
    end
end
end