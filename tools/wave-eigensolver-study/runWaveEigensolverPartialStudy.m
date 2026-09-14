function results = runWaveEigensolverPartialStudy(options)
% Explore a partial generalized solve after exact Dirichlet constraint removal.
arguments
    options.resolutions (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(options.resolutions,4)} = [104 156]
    options.modeCounts (1,:) double {mustBeInteger,mustBePositive} = [8 16 32 64]
    options.kappa (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 3e-4
end
zDomain = [-1000 0];
N2 = @(z) 1e-4*exp(2*z/700);
surface = IMBoundaryCondition(a=0,b=1,c=1,d=0);
rows = cell(0,9);
for nEVP = options.resolutions
    evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=zDomain,k=options.kappa,f0=1e-4, ...
        surfaceBoundary=surface,bottomBoundary=IMBoundaryCondition.dirichlet());
    solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb").configuredForEVP(evp);
    [A,B] = evp.assembleConfigured(solver);
    [scaledA,scaledB,Z,retainedRows,columnScale] = reducedPencil(solver,A,B);
    reducedA = scaledA(retainedRows,:)*Z;
    reducedB = scaledB(retainedRows,:)*Z;
    timer = tic;
    [fullVectors,fullD] = eig(reducedA,reducedB);
    fullSeconds = toc(timer);
    fullLambda = diag(fullD);
    valid = isfinite(fullLambda) & abs(imag(fullLambda)) < 1e-8*max(1,abs(real(fullLambda))) & real(fullLambda) > 0;
    fullLambda = sort(real(fullLambda(valid)),"ascend");
    for nModes = options.modeCounts
        requested = min(nModes,size(reducedA,1)-2);
        eigsOptions = struct("Tolerance",1e-12,"MaxIterations",1000,"Display",false);
        timer = tic;
        try
            [partialVectors,partialD,flag] = eigs(reducedA,reducedB,requested,"smallestabs",eigsOptions);
            partialSeconds = toc(timer);
            partialLambda = diag(partialD);
            valid = isfinite(partialLambda) & abs(imag(partialLambda)) < 1e-8*max(1,abs(real(partialLambda))) & real(partialLambda) > 0;
            partialLambda = sort(real(partialLambda(valid)),"ascend");
            compared = min([requested,numel(partialLambda),numel(fullLambda)]);
            relativeError = max(abs(partialLambda(1:compared)-fullLambda(1:compared))./fullLambda(1:compared));
            reconstructed = columnScale.'.*(Z*partialVectors);
            residual = maximumResidual(A,B,reconstructed,diag(partialD));
            failure = "";
        catch exception
            partialSeconds = toc(timer);
            flag = NaN;
            compared = 0;
            relativeError = NaN;
            residual = NaN;
            failure = string(exception.identifier);
        end
        rows(end+1,:) = {nEVP,requested,compared,flag,relativeError,residual,fullSeconds,partialSeconds,failure}; %#ok<AGROW>
    end
end
results = cell2table(rows,VariableNames=["nEVP","requestedModes","comparedModes","flag", ...
    "maximumRelativeEigenvalueError","maximumPencilResidual","fullSeconds","partialSeconds","failure"]);
end

function [scaledA,scaledB,Z,retainedRows,columnScale] = reducedPencil(solver,A,B)
rowScale = max(max(abs(A),[],2),max(abs(B),[],2));
rowScale(rowScale == 0) = 1;
columnScale = max(1,0:size(A,2)-1).^(-2);
scaledA = (A./rowScale).*columnScale;
scaledB = (B./rowScale).*columnScale;
bottomIndex = solver.boundaryIndex("bottom");
constraint = scaledA(bottomIndex,:);
[~,pivot] = max(abs(constraint));
retainedColumns = setdiff(1:size(A,2),pivot,"stable");
Z = zeros(size(A,2),size(A,2)-1);
Z(retainedColumns,:) = eye(size(A,2)-1);
Z(pivot,:) = -constraint(retainedColumns)/constraint(pivot);
retainedRows = setdiff(1:size(A,1),bottomIndex,"stable");
end

function value = maximumResidual(A,B,vectors,eigenvalues)
values = zeros(1,size(vectors,2));
for iMode = 1:size(vectors,2)
    lambda = eigenvalues(iMode);
    vector = vectors(:,iMode);
    values(iMode) = norm(A*vector-lambda*B*vector,2) ...
        / max((norm(A,2)+abs(lambda)*norm(B,2))*norm(vector,2),realmin);
end
value = max(values);
end
