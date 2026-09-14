function results = measurePartialWaveEigensolver(options)
% Measure full QZ and partial Arnoldi solves on identical wave pencils.
arguments
    options.resolutions (1,:) double {mustBeInteger,mustBePositive} = [104 156]
    options.modeCounts (1,:) double {mustBeInteger,mustBePositive} = [38 64]
    options.kappa (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative} = [0 2*pi/1e5 3e-4 3e-3]
    options.repetitions (1,1) double {mustBeInteger,mustBePositive} = 15
end
profiles = {
    "constant", @(z) 1e-4*ones(size(z));
    "exponential", @(z) 1e-4*exp(2*z/700);
    "strong", @(z) 1e-6+2e-4*exp(z/75)
    };
surface = IMBoundaryCondition(a=0,b=1,c=1,d=0);
rows = cell(0,18);
for iProfile = 1:size(profiles,1)
    for nEVP = options.resolutions
        for kappa = options.kappa
            evp = IMInternalModes.waveModesAtWavenumber(N2=profiles{iProfile,2},zDomain=[-1000 0], ...
                k=kappa,f0=1e-4,surfaceBoundary=surface,bottomBoundary=IMBoundaryCondition.dirichlet());
            solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb").configuredForEVP(evp);
            [A,B] = evp.assembleConfigured(solver);
            for nModes = options.modeCounts
                eigsOptions = deterministicOptions(size(A,1)-1,nModes);
                warmupSolvers(solver,A,B,nModes,eigsOptions);
                fullSeconds = zeros(options.repetitions,4);
                partialSeconds = zeros(options.repetitions,5);
                for repetition = 1:options.repetitions
                    if mod(repetition,2)
                        fullRun = measureFull(A,B);
                        partialRun = measurePartial(solver,A,B,nModes,eigsOptions);
                    else
                        partialRun = measurePartial(solver,A,B,nModes,eigsOptions);
                        fullRun = measureFull(A,B);
                    end
                    fullSeconds(repetition,:) = timingRow(fullRun.timing);
                    partialSeconds(repetition,:) = timingRow(partialRun.timing);
                end
                fullMedian = median(fullSeconds,1);
                partialMedian = median(partialSeconds,1);
                fullLambda = positiveSorted(diag(fullRun.D));
                partialLambda = positiveSorted(diag(partialRun.D));
                compared = min([nModes,numel(fullLambda),numel(partialLambda)]);
                relativeError = max(abs(partialLambda(1:compared)-fullLambda(1:compared))./fullLambda(1:compared));
                residual = maximumResidual(A,B,partialRun.nativeVectors,diag(partialRun.D));
                rows(end+1,:) = {profiles{iProfile,1},nEVP,kappa,nModes,compared,partialRun.flag,relativeError,residual, ...
                    fullMedian(1),fullMedian(2),fullMedian(3),fullMedian(4),partialMedian(1),partialMedian(2), ...
                    partialMedian(3),partialMedian(4),partialMedian(5),partialMedian(5)/fullMedian(4)}; %#ok<AGROW>
            end
        end
    end
end
results = cell2table(rows,VariableNames=["profile","nEVP","kappa","requestedModes","comparedModes", ...
    "flag","maximumRelativeEigenvalueError","maximumPencilResidual","fullPreparationMedianSeconds", ...
    "fullSolveMedianSeconds","fullReconstructionMedianSeconds","fullTotalMedianSeconds", ...
    "partialPreparationMedianSeconds","partialReductionMedianSeconds","partialSolveMedianSeconds", ...
    "partialReconstructionMedianSeconds","partialTotalMedianSeconds","partialTotalOverFullTotal"]);
end

function warmupSolvers(solver,A,B,nModes,eigsOptions)
fullRun = measureFull(A,B);
partialRun = measurePartial(solver,A,B,nModes,eigsOptions);
if size(fullRun.nativeVectors,1) ~= size(A,2) || size(partialRun.nativeVectors,1) ~= size(A,2)
    error("WaveEigensolverStudy:InvalidWarmupOutput", "Warm-up solves did not return native vectors with the expected row count.");
end
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

function options = deterministicOptions(n,nModes)
options = struct("Tolerance",1e-12,"MaxIterations",1000,"Display",false, ...
    "StartVector",ones(n,1)/sqrt(n),"SubspaceDimension",min(n,max(2*nModes+1,20)));
end

function row = timingRow(timing)
if isfield(timing,"reductionSeconds")
    row = [timing.preparationSeconds,timing.reductionSeconds,timing.solveSeconds, ...
        timing.reconstructionSeconds,timing.totalSeconds];
else
    row = [timing.preparationSeconds,timing.solveSeconds,timing.reconstructionSeconds,timing.totalSeconds];
end
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
