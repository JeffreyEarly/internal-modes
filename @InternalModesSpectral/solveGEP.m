% Take matrices A and B from the generalized eigenvalue problem
% (GEP) and returns F,G,h. The last seven arguments are all
% function handles that do as they say.
function [V_cheb,h] = solveGEP(self,A,B,options)

arguments
    self InternalModesSpectral
    A (:,:) double
    B (:,:) double
    options.negativeEigenvalues = 0
    options.minModes = 1
end

if ( any(any(isnan(A))) || any(any(isnan(B))) )
    error('EVP setup fail. Found at least one nan in matrices A and B.\n');
end
[V,D] = eig( A, B );

[h, permutation] = sort(real(self.hFromLambda(diag(D))),'descend');
V_cheb=V(:,permutation);
if options.negativeEigenvalues > 0
    negIndices = find(h<0,options.negativeEigenvalues,'first');
    permutation = cat(1,negIndices,setdiff((1:length(h))',negIndices));
    h = h(permutation);
    V_cheb=V_cheb(:,permutation);
end

resolvedModes = ceil(find(h>0,1,'last')/2); % Have to do ceil, not floor, or we lose the barotropic mode.
if resolvedModes < options.minModes
    ME = MException('GLOceanKit:NeedMorePoints','Returned %d valid modes (%d modes requested) using nEVPs=%d.', resolvedModes, options.minModes, self.nEVP);
    ME = addCause(ME, MException("GLOceanKit:Context", jsonencode(struct("resolvedModes", double(resolvedModes)))));
    throwAsCaller(ME);   % or: throw(ME)
end

end