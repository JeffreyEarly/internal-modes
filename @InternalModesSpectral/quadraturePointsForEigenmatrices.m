function z_g = quadraturePointsForEigenmatrices(self,nPoints,A,B)
% Now we just need to find the roots of the n+1 mode.
% For constant stratification this should give back the
% standard Fourier modes, i.e., an evenly spaced grid.
%
% Note that if the boundary conditions are such that G(0)=0 and
% G(-D)=0, then those two points do not encode any information.
% As such, only the first (nPoints-2) modes will encode any
% useful information. So we'd expect cond(G(:,1:(nPoints-2))))
% to be good (low), but not the next.
if 2*nPoints < self.nEVP
    if ( any(any(isnan(A))) || any(any(isnan(B))) )
        error('GLOceanKit:NaNInMatrix', 'EVP setup fail. Found at least one nan in matrices A and B.\n');
    end
    [V,D] = eig( A, B );

    [h, permutation] = sort(real(self.hFromLambda(diag(D))),'descend');
    G_cheb=V(:,permutation);
    maxModes = ceil(find(h>0,1,'last')/2);

    if maxModes < (nPoints+1)
        error('GLOceanKit:NeedMorePoints', 'Returned %d valid modes (%d quadrature points requested) using nEVPs=%d.',maxModes,nPoints,self.nEVP);
    end

    % Could compute the roots of the F-modes, but nah.
    %                F = self.Diff1_xCheb(G_cheb(:,nPoints-1));
    %                roots = InternalModesSpectral.FindRootsFromChebyshevVector(F(1:end-1), self.z_xLobatto);
    %                z_g = cat(1,min(self.z_xLobatto),reshape(roots,[],1),max(self.z_xLobatto));

    % depending on the boundary conditions and particular
    % problem, the nth mode might contain (n-1), (n), or (n+1)
    % zero crossings. If nPoints are request, we want to include
    % the boundaries in that number.
    if self.upperBoundary == UpperBoundary.mda
        rootMode = nPoints-2;
    elseif self.upperBoundary == UpperBoundary.rigidLid
        % n-th mode has n+1 zeros (including boundaries)
        rootMode = nPoints-1;
    elseif self.upperBoundary == UpperBoundary.freeSurface && self.lowerBoundary == LowerBoundary.noSlip
        % n-th mode has n zeros (including zero at lower
        % boundary, and not zero at upper)
        rootMode = nPoints-1;
    elseif self.upperBoundary == UpperBoundary.freeSurface
        % n-th mode has n zeros (including zero at lower
        % boundary, and not zero at upper)
        rootMode = nPoints;
    end
    rootsVar = InternalModesSpectral.FindRootsFromChebyshevVector(G_cheb(:,rootMode), self.xDomain);

    % First we make sure the roots are within the bounds
    rootsVar(rootsVar<self.xMin) = self.xMin;
    rootsVar(rootsVar>self.xMax) = self.xMax;

    % Add the boundary points---if they are redundant, they will
    % get eliminated below
    rootsVar = cat(1,self.xMin,rootsVar,self.xMax);

    % Then we eliminate any repeats (it happens)
    rootsVar = unique(rootsVar,'stable');

    while (length(rootsVar) > nPoints)
        rootsVar = sort(rootsVar);
        F = InternalModesSpectral.IntegrateChebyshevVector(G_cheb(:,rootMode));
        value = InternalModesSpectral.ValueOfFunctionAtPointOnGrid( rootsVar, self.xDomain, F );
        dv = diff(value);
        [~,minIndex] = min(abs(dv));
        rootsVar(minIndex+1) = [];
    end

    if length(rootsVar) < nPoints
        error('GLOceanKit:NeedMorePoints', 'Returned %d unique roots (requested %d). Maybe need more EVP.', length(rootsVar),nPoints);
    end

    z_g = reshape(rootsVar,[],1);
    z_g = InternalModesSpectral.fInverseBisection(self.x_function,z_g,min(self.zDomain),max(self.zDomain),1e-12);
else
    error('GLOceanKit:NeedMorePoints', 'You need at least twice as many nEVP as points you request');
end
end