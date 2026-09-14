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
rows = cell(0,11);
for iProfile = 1:size(profiles,1)
    for nEVP = options.resolutions
        for kappa = options.kappa
            evp = IMInternalModes.waveModesAtWavenumber(N2=profiles{iProfile,2},zDomain=[-1000 0], ...
                k=kappa,f0=1e-4,surfaceBoundary=surface,bottomBoundary=IMBoundaryCondition.dirichlet());
            solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb").configuredForEVP(evp);
            [A,B] = evp.assembleConfigured(solver);
            [scaledA,scaledB,reducedA,reducedB,reconstruction] = preparePencils(solver,A,B);
            for nModes = options.modeCounts
                eigsOptions = deterministicOptions(size(reducedA,1),nModes);
                eig(scaledA,scaledB);
                eigs(reducedA,reducedB,nModes,"smallestabs",eigsOptions);
                fullSeconds = zeros(options.repetitions,1);
                partialSeconds = zeros(options.repetitions,1);
                for repetition = 1:options.repetitions
                    if mod(repetition,2)
                        fullSeconds(repetition) = measureFull(scaledA,scaledB);
                        partialSeconds(repetition) = measurePartial(reducedA,reducedB,nModes,eigsOptions);
                    else
                        partialSeconds(repetition) = measurePartial(reducedA,reducedB,nModes,eigsOptions);
                        fullSeconds(repetition) = measureFull(scaledA,scaledB);
                    end
                end
                [~,fullD] = eig(scaledA,scaledB);
                [partialVectors,partialD,flag] = eigs(reducedA,reducedB,nModes,"smallestabs",eigsOptions);
                fullLambda = positiveSorted(diag(fullD));
                partialLambda = positiveSorted(diag(partialD));
                compared = min([nModes,numel(fullLambda),numel(partialLambda)]);
                relativeError = max(abs(partialLambda(1:compared)-fullLambda(1:compared))./fullLambda(1:compared));
                nativeVectors = reconstruction*partialVectors;
                residual = maximumResidual(A,B,nativeVectors,diag(partialD));
                rows(end+1,:) = {profiles{iProfile,1},nEVP,kappa,nModes,compared,flag,relativeError,residual, ...
                    median(fullSeconds),median(partialSeconds),median(partialSeconds)/median(fullSeconds)}; %#ok<AGROW>
            end
        end
    end
end
results = cell2table(rows,VariableNames=["profile","nEVP","kappa","requestedModes","comparedModes", ...
    "flag","maximumRelativeEigenvalueError","maximumPencilResidual","fullMedianSeconds", ...
    "partialMedianSeconds","partialOverFull"]);
end

function [scaledA,scaledB,reducedA,reducedB,reconstruction] = preparePencils(solver,A,B)
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
reducedA = scaledA(retainedRows,:)*Z;
reducedB = scaledB(retainedRows,:)*Z;
reconstruction = columnScale.'.*Z;
end

function options = deterministicOptions(n,nModes)
options = struct("Tolerance",1e-12,"MaxIterations",1000,"Display",false, ...
    "StartVector",ones(n,1)/sqrt(n),"SubspaceDimension",min(n,max(2*nModes+1,20)));
end

function seconds = measureFull(A,B)
timer = tic;
eig(A,B);
seconds = toc(timer);
end

function seconds = measurePartial(A,B,nModes,options)
timer = tic;
eigs(A,B,nModes,"smallestabs",options);
seconds = toc(timer);
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
