function results = runWaveEigensolverCorrectnessStudy(options)
% Compare the production generalized solve with an inverse-pencil prototype.
arguments
    options.resolutions (1,:) double {mustBeInteger,mustBeGreaterThanOrEqual(options.resolutions,4)} = [69 104 156]
    options.kappa (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative} = [0 2*pi/1e5 3e-4 3e-3]
    options.nModes (1,1) double {mustBeInteger,mustBePositive} = 64
    options.coordinateKind (1,1) string {mustBeMember(options.coordinateKind,["z","wkb"])} = "wkb"
end

zDomain = [-1000 0];
f0 = 1e-4;
g = 9.81;
surface = IMBoundaryCondition(a=0,b=1,c=1,d=0);
bottom = IMBoundaryCondition.dirichlet();
profiles = {
    "constant", @(z) 1e-4*ones(size(z));
    "exponential", @(z) 1e-4*exp(2*z/700);
    "strong", @(z) 1e-6+2e-4*exp(z/75)
    };

rows = cell(0,16);
for iProfile = 1:size(profiles,1)
    profileName = profiles{iProfile,1};
    N2 = profiles{iProfile,2};
    for nEVP = options.resolutions
        for kappa = options.kappa
            evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=zDomain,k=kappa,f0=f0,g=g, ...
                surfaceBoundary=surface,bottomBoundary=bottom);
            solver = IMSolverSpectral(nEVP=nEVP,coordinateKind=options.coordinateKind).configuredForEVP(evp);
            [A,B,samples] = evp.assembleConfigured(solver);
            diagnostics = evp.preparedModeSelectionDiagnostics(samples,A);
            baseline = solver.solveEVP(evp,nModes=options.nModes);
            candidates = {
                "inverseA", @() inversePencilSolve(solver,evp,A,B,diagnostics,options.nModes);
                "reducedMetric", @() reducedMetricSolve(solver,evp,A,B,diagnostics,options.nModes);
                "gridMetric", @() gridMetricSolve(solver,evp,A,B,diagnostics,options.nModes)
                };
            for iCandidate = 1:size(candidates,1)
                [candidate,metrics] = candidates{iCandidate,2}();
                z = linspace(zDomain(1),zDomain(2),257).';
                baselineG = baseline.G(z);
                candidateG = candidate.G(z);
                baselineF = baseline.F(z);
                candidateF = candidate.F(z);
                signs = sign(sum(baselineG.*candidateG,1));
                signs(signs == 0) = 1;
                candidateG = candidateG.*signs;
                candidateF = candidateF.*signs;
                relativeH = max(abs(candidate.h-baseline.h)./max(abs(baseline.h),realmin));
                relativeG = norm(candidateG-baselineG,"fro")/max(norm(baselineG,"fro"),realmin);
                relativeF = norm(candidateF-baselineF,"fro")/max(norm(baselineF,"fro"),realmin);
                labelsEqual = isequal(candidate.modeNumber,baseline.modeNumber);
                rows(end+1,:) = {profileName,candidates{iCandidate,1},nEVP,kappa,metrics.conditionIndicator, ...
                    metrics.metricRank,metrics.finiteCandidateCount,metrics.maximumPencilResidual, ...
                    metrics.minimumAbsoluteEigenvalue,relativeH,relativeG,relativeF,labelsEqual, ...
                    metrics.generalizedSeconds,metrics.candidateSeconds,metrics.candidateOverGeneralized}; %#ok<AGROW>
            end
        end
    end
end
results = cell2table(rows,VariableNames=["profile","candidate","nEVP","kappa","conditionIndicator","metricRank", ...
    "finiteCandidateCount","maximumPencilResidual","minimumAbsoluteEigenvalue","maximumRelativeHError", ...
    "relativeGError","relativeFError","labelsEqual","generalizedSeconds","candidateSeconds", ...
    "candidateOverGeneralized"]);
end

function [basis,metrics] = inversePencilSolve(solver,evp,A,B,diagnostics,nModes)
rowScale = max(max(abs(A),[],2),max(abs(B),[],2));
rowScale(rowScale == 0) = 1;
solveA = A./rowScale;
solveB = B./rowScale;
columnScale = max(1,0:size(A,2)-1).^(-2);
scaledA = solveA.*columnScale;
scaledB = solveB.*columnScale;

timer = tic;
[generalizedVectors,generalizedValues] = eig(scaledA,scaledB); %#ok<ASGLU>
generalizedSeconds = toc(timer);
timer = tic;
factorization = decomposition(scaledA);
inversePencil = factorization\scaledB;
[vectors,H] = eig(inversePencil);
inverseSeconds = toc(timer);
h = diag(H);
lambda = 1./h;
vectors = columnScale.'.*vectors;
valid = isfinite(real(lambda)) & isfinite(imag(lambda)) ...
    & abs(imag(lambda)) < 1e-8*max(1,abs(real(lambda)));
vectors = real(vectors(:,valid));
lambda = real(lambda(valid));
familyValid = evp.finiteGeneralizedEigenpairMask(vectors,solveB);
vectors = vectors(:,familyValid);
lambda = lambda(familyValid);
selection = evp.selectModes(lambda(:),nModes,solver,A,diagnostics=diagnostics);
lambda = lambda(selection.sortIndex);
lambda(selection.modeNumber == 0) = 0;
vectors = vectors(:,selection.sortIndex);
basis = evp.makeBasisSet(solver,vectors,lambda(:).',selection.modeNumber,selection.modeSelectionDiagnostics);
basis = basis.orientModeSigns();

residuals = zeros(1,size(vectors,2));
normA = norm(A,2);
normB = norm(B,2);
for iMode = 1:size(vectors,2)
    vector = vectors(:,iMode);
    denominator = (normA+abs(lambda(iMode))*normB)*norm(vector,2);
    residuals(iMode) = norm(A*vector-lambda(iMode)*B*vector,2)/max(denominator,realmin);
end
finiteH = abs(h(isfinite(h) & h ~= 0));
metrics = struct("conditionIndicator",rcond(scaledA),"metricRank",rank(scaledB), ...
    "finiteCandidateCount",numel(lambda),"maximumPencilResidual",max(residuals), ...
    "minimumAbsoluteEigenvalue",min(finiteH),"generalizedSeconds",generalizedSeconds, ...
    "candidateSeconds",inverseSeconds,"candidateOverGeneralized",inverseSeconds/generalizedSeconds);
end

function [basis,metrics] = reducedMetricSolve(solver,evp,A,B,diagnostics,nModes)
rowScale = max(max(abs(A),[],2),max(abs(B),[],2));
rowScale(rowScale == 0) = 1;
solveA = A./rowScale;
solveB = B./rowScale;
columnScale = max(1,0:size(A,2)-1).^(-2);
scaledA = solveA.*columnScale;
scaledB = solveB.*columnScale;
timer = tic;
eig(scaledA,scaledB);
generalizedSeconds = toc(timer);

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
timer = tic;
metricFactorization = decomposition(reducedB);
standardMatrix = metricFactorization\reducedA;
[reducedVectors,D] = eig(standardMatrix);
candidateSeconds = toc(timer);
lambda = diag(D);
vectors = columnScale.'.*(Z*reducedVectors);
valid = isfinite(real(lambda)) & isfinite(imag(lambda)) ...
    & abs(imag(lambda)) < 1e-8*max(1,abs(real(lambda)));
vectors = real(vectors(:,valid));
lambda = real(lambda(valid));
familyValid = evp.finiteGeneralizedEigenpairMask(vectors,solveB);
vectors = vectors(:,familyValid);
lambda = lambda(familyValid);
selection = evp.selectModes(lambda(:),nModes,solver,A,diagnostics=diagnostics);
lambda = lambda(selection.sortIndex);
lambda(selection.modeNumber == 0) = 0;
vectors = vectors(:,selection.sortIndex);
basis = evp.makeBasisSet(solver,vectors,lambda(:).',selection.modeNumber,selection.modeSelectionDiagnostics);
basis = basis.orientModeSigns();

residuals = zeros(1,size(vectors,2));
normA = norm(A,2);
normB = norm(B,2);
for iMode = 1:size(vectors,2)
    vector = vectors(:,iMode);
    denominator = (normA+abs(lambda(iMode))*normB)*norm(vector,2);
    residuals(iMode) = norm(A*vector-lambda(iMode)*B*vector,2)/max(denominator,realmin);
end
finiteLambda = abs(lambda(isfinite(lambda) & lambda ~= 0));
metrics = struct("conditionIndicator",rcond(reducedB),"metricRank",rank(reducedB), ...
    "finiteCandidateCount",numel(lambda),"maximumPencilResidual",max(residuals), ...
    "minimumAbsoluteEigenvalue",min(finiteLambda),"generalizedSeconds",generalizedSeconds, ...
    "candidateSeconds",candidateSeconds,"candidateOverGeneralized",candidateSeconds/generalizedSeconds);
end

function [basis,metrics] = gridMetricSolve(solver,evp,A,B,diagnostics,nModes)
rowScale = max(max(abs(A),[],2),max(abs(B),[],2));
rowScale(rowScale == 0) = 1;
columnScale = max(1,0:size(A,2)-1).^(-2);
timer = tic;
eig((A./rowScale).*columnScale,(B./rowScale).*columnScale);
generalizedSeconds = toc(timer);

bottomIndex = solver.boundaryIndex("bottom");
retained = setdiff(1:size(A,1),bottomIndex,"stable");
timer = tic;
gridA = A/solver.T;
gridB = B/solver.T;
reducedA = gridA(retained,retained);
reducedB = gridB(retained,retained);
[gridVectors,D] = eig(reducedB\reducedA);
candidateSeconds = toc(timer);
lambda = diag(D);
gridModes = zeros(size(A,1),size(gridVectors,2));
gridModes(retained,:) = gridVectors;
vectors = solver.T\gridModes;
valid = isfinite(real(lambda)) & isfinite(imag(lambda)) ...
    & abs(imag(lambda)) < 1e-8*max(1,abs(real(lambda)));
vectors = real(vectors(:,valid));
lambda = real(lambda(valid));
familyValid = evp.finiteGeneralizedEigenpairMask(vectors,B./rowScale);
vectors = vectors(:,familyValid);
lambda = lambda(familyValid);
selection = evp.selectModes(lambda(:),nModes,solver,A,diagnostics=diagnostics);
lambda = lambda(selection.sortIndex);
lambda(selection.modeNumber == 0) = 0;
vectors = vectors(:,selection.sortIndex);
basis = evp.makeBasisSet(solver,vectors,lambda(:).',selection.modeNumber,selection.modeSelectionDiagnostics);
basis = basis.orientModeSigns();

residuals = zeros(1,size(vectors,2));
normA = norm(A,2);
normB = norm(B,2);
for iMode = 1:size(vectors,2)
    vector = vectors(:,iMode);
    denominator = (normA+abs(lambda(iMode))*normB)*norm(vector,2);
    residuals(iMode) = norm(A*vector-lambda(iMode)*B*vector,2)/max(denominator,realmin);
end
finiteLambda = abs(lambda(isfinite(lambda) & lambda ~= 0));
metrics = struct("conditionIndicator",rcond(reducedB),"metricRank",rank(reducedB), ...
    "finiteCandidateCount",numel(lambda),"maximumPencilResidual",max(residuals), ...
    "minimumAbsoluteEigenvalue",min(finiteLambda),"generalizedSeconds",generalizedSeconds, ...
    "candidateSeconds",candidateSeconds,"candidateOverGeneralized",candidateSeconds/generalizedSeconds);
end
