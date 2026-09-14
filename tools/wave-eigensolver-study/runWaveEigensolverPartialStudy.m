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
rows = cell(0,17);
for nEVP = options.resolutions
    evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=zDomain,k=options.kappa,f0=1e-4, ...
        surfaceBoundary=surface,bottomBoundary=IMBoundaryCondition.dirichlet());
    solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb").configuredForEVP(evp);
    [A,B] = evp.assembleConfigured(solver);
    fullRun = measureFull(A,B);
    fullLambda = positiveSorted(diag(fullRun.D));
    for nModes = options.modeCounts
        requested = min(nModes,size(A,1)-3);
        eigsOptions = struct("Tolerance",1e-12,"MaxIterations",1000,"Display",false);
        partialTimer = tic;
        try
            partialRun = measurePartial(solver,A,B,requested,eigsOptions);
            partialLambda = positiveSorted(diag(partialRun.D));
            compared = min([requested,numel(partialLambda),numel(fullLambda)]);
            relativeError = max(abs(partialLambda(1:compared)-fullLambda(1:compared))./fullLambda(1:compared));
            residual = maximumResidual(A,B,partialRun.nativeVectors,diag(partialRun.D));
            flag = partialRun.flag;
            partialTiming = partialRun.timing;
            failure = "";
        catch exception
            partialTiming = struct("preparationSeconds",NaN,"reductionSeconds",NaN,"solveSeconds",NaN, ...
                "reconstructionSeconds",NaN,"totalSeconds",toc(partialTimer));
            flag = NaN;
            compared = 0;
            relativeError = NaN;
            residual = NaN;
            failure = string(exception.identifier);
        end
        rows(end+1,:) = {nEVP,requested,compared,flag,relativeError,residual, ...
            fullRun.timing.preparationSeconds,fullRun.timing.solveSeconds, ...
            fullRun.timing.reconstructionSeconds,fullRun.timing.totalSeconds, ...
            partialTiming.preparationSeconds,partialTiming.reductionSeconds,partialTiming.solveSeconds, ...
            partialTiming.reconstructionSeconds,partialTiming.totalSeconds, ...
            partialTiming.totalSeconds/fullRun.timing.totalSeconds,failure}; %#ok<AGROW>
    end
end
results = cell2table(rows,VariableNames=["nEVP","requestedModes","comparedModes","flag", ...
    "maximumRelativeEigenvalueError","maximumPencilResidual","fullPreparationSeconds", ...
    "fullSolveSeconds","fullReconstructionSeconds","fullTotalSeconds","partialPreparationSeconds", ...
    "partialReductionSeconds","partialSolveSeconds","partialReconstructionSeconds", ...
    "partialTotalSeconds","partialTotalOverFullTotal","failure"]);
end

function run = measureFull(A,B)
preparationTimer = tic;
[scaledA,scaledB,columnScale] = scalePencil(A,B);
preparationSeconds = toc(preparationTimer);
solveTimer = tic;
[scaledVectors,D] = eig(scaledA,scaledB);
solveSeconds = toc(solveTimer);
reconstructionTimer = tic;
nativeVectors = columnScale.'.*scaledVectors;
reconstructionSeconds = toc(reconstructionTimer);
run = struct("nativeVectors",nativeVectors,"D",D,"timing",struct( ...
    "preparationSeconds",preparationSeconds,"solveSeconds",solveSeconds, ...
    "reconstructionSeconds",reconstructionSeconds, ...
    "totalSeconds",preparationSeconds+solveSeconds+reconstructionSeconds));
end

function run = measurePartial(solver,A,B,nModes,eigsOptions)
preparationTimer = tic;
[scaledA,scaledB,columnScale] = scalePencil(A,B);
preparationSeconds = toc(preparationTimer);
reductionTimer = tic;
[reducedA,reducedB,reconstruction] = reducePencil(solver,scaledA,scaledB,columnScale);
reductionSeconds = toc(reductionTimer);
solveTimer = tic;
[reducedVectors,D,flag] = eigs(reducedA,reducedB,nModes,"smallestabs",eigsOptions);
solveSeconds = toc(solveTimer);
reconstructionTimer = tic;
nativeVectors = reconstruction*reducedVectors;
reconstructionSeconds = toc(reconstructionTimer);
run = struct("nativeVectors",nativeVectors,"D",D,"flag",flag,"timing",struct( ...
    "preparationSeconds",preparationSeconds,"reductionSeconds",reductionSeconds, ...
    "solveSeconds",solveSeconds,"reconstructionSeconds",reconstructionSeconds, ...
    "totalSeconds",preparationSeconds+reductionSeconds+solveSeconds+reconstructionSeconds));
end

function [scaledA,scaledB,columnScale] = scalePencil(A,B)
rowScale = max(max(abs(A),[],2),max(abs(B),[],2));
rowScale(rowScale == 0) = 1;
columnScale = max(1,0:size(A,2)-1).^(-2);
scaledA = (A./rowScale).*columnScale;
scaledB = (B./rowScale).*columnScale;
end

function [reducedA,reducedB,reconstruction] = reducePencil(solver,scaledA,scaledB,columnScale)
bottomIndex = solver.boundaryIndex("bottom");
constraint = scaledA(bottomIndex,:);
[~,pivot] = max(abs(constraint));
retainedColumns = setdiff(1:size(scaledA,2),pivot,"stable");
Z = zeros(size(scaledA,2),size(scaledA,2)-1);
Z(retainedColumns,:) = eye(size(scaledA,2)-1);
Z(pivot,:) = -constraint(retainedColumns)/constraint(pivot);
retainedRows = setdiff(1:size(scaledA,1),bottomIndex,"stable");
reducedA = scaledA(retainedRows,:)*Z;
reducedB = scaledB(retainedRows,:)*Z;
reconstruction = columnScale.'.*Z;
end

function values = positiveSorted(values)
valid = isfinite(values) & abs(imag(values)) < 1e-8*max(1,abs(real(values))) & real(values) > 0;
values = sort(real(values(valid)),"ascend");
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
