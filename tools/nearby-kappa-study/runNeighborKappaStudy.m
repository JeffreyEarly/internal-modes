function summary = runNeighborKappaStudy(outputDirectory,options)
% Offline #21 research. No production reuse policy or defaults are changed.
arguments
    outputDirectory (1,1) string
    options.caseIds (1,:) double {mustBeInteger,mustBePositive} = 1:96
    options.maxSeconds (1,1) double {mustBePositive} = 600
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
configuration = studyConfiguration();
cases = configuration.cases;
if any(options.caseIds > height(cases)), error('IM21:InvalidCase','caseIds must select declared cases.'); end
cases = cases(options.caseIds,:);
configuration.selectedCaseIds = cases.caseId;
configuration.maxSeconds = options.maxSeconds;
save(fullfile(outputDirectory,'configuration.mat'),'configuration');
records = struct([]);
runTimer = tic;
for environment = unique(cases.environment,'stable').'
    selected = cases(cases.environment == environment,:);
    problem = studyProblem(selected.profile(1),selected.surface(1),configuration);
    requestedKappa = reshape(unique([selected.sourceKappa;selected.requestedKappa],'stable'),1,[]);
    baselines = cell(1,3); constructionCosts = cell(1,3);
    for iResolution = 1:3
        solver = IMSolverSpectral(nEVP=configuration.orders(iResolution),coordinateKind='wkb');
        [baselines{iResolution},constructionCosts{iResolution}] = solver.solveWaveModesAtWavenumbers(requestedKappa,N2=problem.N2,zDomain=problem.zDomain,f0=problem.f0,g=problem.g,surfaceBoundary=problem.surfaceBoundary,nModes=configuration.referenceModes);
    end
    retained = struct('cases',selected,'problem',problem,'resolutions',configuration.orders,'collections',{baselines},'constructionCosts',{constructionCosts},'candidateRecords',{{}});
    for iCase = 1:height(selected)
        if toc(runTimer) > options.maxSeconds
            retained.stoppedForBudget = true;
            save(fullfile(outputDirectory,sprintf('environment-%02d.mat',environment)),'retained','-v7.3');
            summary = finish(outputDirectory,records,configuration,toc(runTimer),false);
            return
        end
        row = selected(iCase,:);
        targetIndex = find(requestedKappa == row.requestedKappa,1);
        sourceIndex = find(requestedKappa == row.sourceKappa,1);
        exact = baselines{1}.bases{targetIndex};
        neighbor = baselines{1}.bases{sourceIndex};
        coarse = baselines{2}.bases{targetIndex};
        reference = baselines{3}.bases{targetIndex};
        referenceTimer = tic;
        comparison = compareBases(coarse,reference,row.requestedKappa,configuration);
        residual = physicalResidual(reference,row.requestedKappa,configuration);
        referenceQualified = referencePasses(comparison,residual,configuration);
        escalationSeconds = 0;
        if ~referenceQualified
            timer = tic;
            escalated = IMSolverSpectral(nEVP=configuration.escalationOrder,coordinateKind='wkb').solveEVP(reference.evp,nModes=configuration.referenceModes);
            escalationSeconds = toc(timer);
            comparison = compareBases(reference,escalated,row.requestedKappa,configuration);
            reference = escalated;
            residual = physicalResidual(reference,row.requestedKappa,configuration);
            referenceQualified = referencePasses(comparison,residual,configuration);
        end
        referenceValidationSeconds = toc(referenceTimer)-escalationSeconds;
        analytic = analyticalControl(problem,reference,row.requestedKappa,configuration);
        referenceValidationSeconds = referenceValidationSeconds+analytic.seconds;
        if analytic.available
            referenceQualified = referenceQualified && referencePasses(analytic.comparison,residual,configuration);
        end
        baselineTimer = tic;
        baselineComparison = compareBases(exact,reference,row.requestedKappa,configuration);
        baselineResidual = physicalResidual(exact,row.requestedKappa,configuration);
        baselineValidationSeconds = toc(baselineTimer);
        setupTimer = tic;
        solver = exact.solver;
        [A,B,samples] = exact.evp.assembleConfigured(solver);
        requestedAssemblySeconds = toc(setupTimer);
        setupTimer = tic;
        selectionDiagnostics = exact.evp.preparedModeSelectionDiagnostics(samples,A);
        requestedSelectionSeconds = toc(setupTimer);
        [warm,warmCost] = acceleratedSolve(exact.evp,solver,A,B,neighbor.nativeModes,selectionDiagnostics,configuration,true);
        [cold,coldCost] = acceleratedSolve(exact.evp,solver,A,B,neighbor.nativeModes,selectionDiagnostics,configuration,false);
        validationTimer = tic;
        if isempty(warm)
            warmComparison = failedComparison(configuration);
            warmResidual = failedResidual();
            warmToFull = failedComparison(configuration);
        else
            warmComparison = compareBases(warm,reference,row.requestedKappa,configuration);
            warmResidual = physicalResidual(warm,row.requestedKappa,configuration);
            warmToFull = compareBases(warm,exact,row.requestedKappa,configuration);
        end
        warmValidationSeconds = toc(validationTimer);
        if isempty(cold)
            coldToFull = failedComparison(configuration);
        else
            coldToFull = compareBases(cold,exact,row.requestedKappa,configuration);
        end
        approximationTimer = tic;
        approximateComparison = compareBases(neighbor,reference,row.requestedKappa,configuration);
        approximateResidual = physicalResidual(neighbor,row.requestedKappa,configuration);
        approximationValidationSeconds = toc(approximationTimer);
        fullCost = constructionCosts{1}.solves(targetIndex,:);
        fullSharedSeconds = (constructionCosts{1}.setupSeconds+constructionCosts{1}.sharedAssemblySeconds+constructionCosts{1}.collectionSeconds)/numel(requestedKappa);
        fullSolveSeconds = fullCost.assemblySeconds+fullCost.eigensolveSeconds+fullCost.finalizationSeconds+fullSharedSeconds;
        sourceCost = constructionCosts{1}.solves(sourceIndex,:);
        neighborAcquisitionSeconds = sourceCost.assemblySeconds+sourceCost.eigensolveSeconds+sourceCost.finalizationSeconds+fullSharedSeconds;
        warmMatchesFull = warmCost.flag == 0 && all(warmToFull.eigenvalueError <= configuration.candidateEigenvalueTolerance) && all(warmToFull.subspaceError <= configuration.candidateSubspaceTolerance);
        fallbackRequired = ~warmMatchesFull || ~referenceQualified || ~referencePasses(warmComparison,warmResidual,configuration);
        % A full solve is already retained as the independent comparison.
        % Charge its measured cost on fallback; never claim it repairs an
        % under-resolved full-solve reference or missing product validation.
        warmTotalWithFallbackSeconds = requestedAssemblySeconds+requestedSelectionSeconds+warmCost.totalSeconds+warmValidationSeconds+fallbackRequired*fullSolveSeconds;
        referenceSolveSeconds = 0;
        for iResolution = 2:3
            cost = constructionCosts{iResolution};
            referenceSolveSeconds = referenceSolveSeconds+cost.solves.eigensolveSeconds(targetIndex)+cost.solves.finalizationSeconds(targetIndex)+cost.solves.assemblySeconds(targetIndex)+(cost.setupSeconds+cost.sharedAssemblySeconds+cost.collectionSeconds)/numel(requestedKappa);
        end
        record = struct('caseId',row.caseId,'profile',row.profile,'surface',row.surface,'sourceKappa',row.sourceKappa,'requestedKappa',row.requestedKappa,'relativeOffset',row.relativeOffset, ...
            'referenceQualifiedExperimental',referenceQualified,'referenceOrder',reference.solver.nEVP,'referenceComparison',comparison,'referenceResidual',residual, ...
            'baselineComparison',baselineComparison,'baselineResidual',baselineResidual,'warmComparison',warmComparison,'warmResidual',warmResidual,'warmToFull',warmToFull,'coldToFull',coldToFull, ...
            'approximateComparison',approximateComparison,'approximateResidual',approximateResidual,'warmCost',warmCost,'coldCost',coldCost,'analyticControl',analytic, ...
            'fullSolveSeconds',fullSolveSeconds,'fullSharedPreparationAmortizedSeconds',fullSharedSeconds,'neighborAcquisitionSeconds',neighborAcquisitionSeconds,'neighborStateAssumption','existing solved neighbor; acquisition reported separately','requestedAssemblySeconds',requestedAssemblySeconds,'requestedSelectionSeconds',requestedSelectionSeconds,'warmValidationSeconds',warmValidationSeconds,'baselineValidationSeconds',baselineValidationSeconds, ...
            'referenceSolveSeconds',referenceSolveSeconds+escalationSeconds,'referenceValidationSeconds',referenceValidationSeconds,'approximationValidationSeconds',approximationValidationSeconds, ...
            'warmMatchesFull',warmMatchesFull,'fallbackRequired',fallbackRequired,'warmTotalWithFallbackSeconds',warmTotalWithFallbackSeconds, ...
            'validatedWithReferenceSeconds',warmTotalWithFallbackSeconds+referenceSolveSeconds+escalationSeconds+referenceValidationSeconds+(~fallbackRequired)*fullSolveSeconds+baselineValidationSeconds, ...
            'validationUsesStoredOracles',true, ...
            'productAssessmentStatus','deferred pending #19','productionQualification',false);
        if isempty(records)
            records = record;
        else
            records(end+1) = record; %#ok<AGROW>
        end
        retained.candidateRecords{end+1} = struct('caseId',row.caseId,'sourceBasis',neighbor,'exactBasis',exact,'referenceBasis',reference,'warmBasis',warm,'coldBasis',cold,'record',record);
        fprintf('case %d profile=%s surface=%s delta=%g ref=%d qualified=%d warmFlag=%d fallback=%d approxH1=%g elapsed=%.1fs\n',row.caseId,row.profile,row.surface,row.relativeOffset,reference.solver.nEVP,referenceQualified,warmCost.flag,fallbackRequired,max(approximateComparison.subspaceError),toc(runTimer));
        save(fullfile(outputDirectory,'progress.mat'),'records','configuration');
    end
    component = whos('retained');
    retained.workspaceBytesEstimate = component.bytes;
    save(fullfile(outputDirectory,sprintf('environment-%02d.mat',environment)),'retained','-v7.3');
end
summary = finish(outputDirectory,records,configuration,toc(runTimer),true);
end

function c = studyConfiguration()
c = struct('orders',[64 96 144],'escalationOrder',216,'bands',[8 16],'referenceModes',20,'D',1000,'f0',1e-4,'g',9.81, ...
    'referenceEigenvalueTolerance',1e-7,'referenceSubspaceTolerance',1e-5,'equationTolerance',1e-6,'boundaryTolerance',1e-8, ...
    'candidateEigenvalueTolerance',1e-8,'candidateSubspaceTolerance',1e-6,'eigsTolerance',1e-10,'maxIterations',120,'krylovDimension',48,'referenceGridCount',513);
c.guardMeaning = 'Experimental comparisons only; not accepted #10 diagnostics or production acceptance guarantees.';
c.missingQualification = 'Independent pressure/momentum residual, source products and sampled transforms are not qualified.';
profiles = ["constant","exponential","pycnocline","doubleSharp"];
surfaces = ["rigid","free"];
anchors = [0.1 1 10]/c.D;
offsets = [0 1e-3 1e-2 1e-1];
nCases = numel(profiles)*numel(surfaces)*numel(anchors)*numel(offsets);
caseId = (1:nCases).'; environment = zeros(nCases,1); profile = strings(nCases,1); surface = strings(nCases,1); sourceKappa = zeros(nCases,1); requestedKappa = zeros(nCases,1); relativeOffset = zeros(nCases,1);
iCase = 0;
for iProfile = 1:4
    for iSurface = 1:2
        for anchor = anchors
            for offset = offsets
                iCase = iCase+1;
                environment(iCase,1) = 2*(iProfile-1)+iSurface;
                profile(iCase,1) = profiles(iProfile);
                surface(iCase,1) = surfaces(iSurface);
                sourceKappa(iCase,1) = anchor;
                requestedKappa(iCase,1) = anchor*(1+offset);
                relativeOffset(iCase,1) = offset;
            end
        end
    end
end
c.cases = table(caseId,environment,profile,surface,sourceKappa,requestedKappa,relativeOffset);
end

function problem = studyProblem(profile,surface,c)
D = c.D;
switch profile
    case 'constant'
        N2 = @(z) 1e-4*ones(size(z));
    case 'exponential'
        N2 = @(z) 1e-4*exp(2*z/D);
    case 'pycnocline'
        N2 = @(z) 1e-5+9e-5*exp(-((z+0.35*D)/(0.08*D)).^2);
    case 'doubleSharp'
        N2 = @(z) 1e-5+9e-5*(exp(-((z+0.3*D)/(0.02*D)).^2)+exp(-((z+0.7*D)/(0.02*D)).^2));
end
boundary = IMBoundaryCondition.dirichlet();
if surface == "free", boundary = IMBoundaryCondition(a=0,b=1,c=1,d=0); end
problem = struct('profile',profile,'surface',surface,'N2',N2,'zDomain',[-D 0],'f0',c.f0,'g',c.g,'surfaceBoundary',boundary);
end

function [basis,cost] = acceleratedSolve(evp,solver,A,B,neighborNative,selectionDiagnostics,c,warm)
totalTimer = tic;
factorTimer = tic;
factor = decomposition(A,'lu');
cost = struct('factorizationSeconds',toc(factorTimer),'eigensolveSeconds',0,'finalizationSeconds',0,'operatorApplications',0,'flag',-1,'failure',"",'totalSeconds',0);
applicationCount = 0;
n = size(A,1);
if warm
    v0 = neighborNative*(1./(1:size(neighborNative,2))).';
else
    v0 = sin((1:n).'*sqrt(2))+cos((1:n).'*sqrt(3));
end
v0 = v0/norm(v0);
basis = [];
try
    timer = tic;
    options = struct('tol',c.eigsTolerance,'maxit',c.maxIterations,'p',min(n,c.krylovDimension),'v0',v0,'isreal',true,'issym',false,'disp',0);
    [V,D,flag] = eigs(@apply,n,c.referenceModes,'largestabs',options);
    cost.eigensolveSeconds = toc(timer);
    cost.flag = flag;
    timer = tic;
    lambda = 1./diag(D);
    valid = isfinite(real(lambda)) & isfinite(imag(lambda)) & abs(imag(lambda)) < 1e-8*max(1,abs(real(lambda)));
    V = real(V(:,valid)); lambda = real(lambda(valid));
    valid = evp.finiteGeneralizedEigenpairMask(V,B);
    V = V(:,valid); lambda = lambda(valid);
    selection = evp.selectModes(lambda,c.referenceModes,solver,A,diagnostics=selectionDiagnostics);
    if numel(selection.modeNumber) < c.referenceModes
        error('IM21:InsufficientRitzModes','Insufficient finite real Ritz modes.');
    end
    eigenvalues = lambda(selection.sortIndex);
    eigenvalues(selection.modeNumber == 0) = 0;
    basis = evp.makeBasisSet(solver,V(:,selection.sortIndex),eigenvalues(:).',selection.modeNumber,selection.modeSelectionDiagnostics);
    basis = basis.orientModeSigns();
    cost.finalizationSeconds = toc(timer);
catch exception
    cost.flag = -1;
    cost.failure = string(exception.identifier)+": "+string(exception.message);
end
cost.operatorApplications = applicationCount;
cost.totalSeconds = toc(totalTimer);
    function y = apply(x)
        applicationCount = applicationCount+size(x,2);
        y = factor\(B*x);
    end
end

function result = compareBases(candidate,reference,requestedKappa,c)
z = linspace(-c.D,0,c.referenceGridCount).';
w = [0.5;ones(numel(z)-2,1);0.5]*(c.D/(numel(z)-1));
[G,Gz,F] = fields(candidate,z);
[Gr,Gzr,Fr] = fields(reference,z);
result = failedComparison(c);
for iBand = 1:numel(c.bands)
    m = c.bands(iBand);
    lambda = candidate.eigenvalues(1:m); lambdaReference = reference.eigenvalues(1:m);
    result.cutoffRelativeGap(iBand) = abs(reference.eigenvalues(m+1)-reference.eigenvalues(m))/max(abs(reference.eigenvalues(m:m+1)));
    result.eigenvalueError(iBand) = max(abs(lambda-lambdaReference)./max(abs(lambdaReference),eps(max(abs(lambdaReference)))));
    omega = sqrt(c.g*requestedKappa^2*candidate.h(1:m)+c.f0^2);
    omegaReference = sqrt(c.g*requestedKappa^2*reference.h(1:m)+c.f0^2);
    result.frequencyError(iBand) = max(abs(omega-omegaReference)./max(abs(omegaReference),c.f0));
    left = [sqrt(w/c.D).*G(:,1:m);sqrt(w*c.D).*Gz(:,1:m)];
    right = [sqrt(w/c.D).*Gr(:,1:m);sqrt(w*c.D).*Gzr(:,1:m)];
    [Q,~] = qr(left,0); [Qr,~] = qr(right,0);
    % Residual projector form avoids cancellation in sqrt(1-sigma_min^2).
    result.subspaceError(iBand) = norm(Q-Qr*(Qr'*Q),2);
    signs = sign(sum(left.*right,1)); signs(signs == 0) = 1;
    result.individualH1Error(iBand) = max(vecnorm(left.*signs-right)./max(vecnorm(right),realmin));
    result.FError(iBand) = max(vecnorm(sqrt(w).*F(:,1:m).*signs-sqrt(w).*Fr(:,1:m))./max(vecnorm(sqrt(w).*Fr(:,1:m)),realmin));
end
lambda = reference.eigenvalues;
result.minimumRelativeGap = min(abs(diff(lambda(1:c.bands(end))))./max(abs(lambda(1:c.bands(end)-1)),abs(lambda(2:c.bands(end)))));
result.candidateLabels = candidate.modeNumber;
result.referenceLabels = reference.modeNumber;
result.matching = 'Same numerical ordering by eigenvalue; ordinal labels retained. Cluster ambiguity is not resolved by column errors.';
end

function [G,Gz,F] = fields(basis,z)
factors = basis.normalizationFactors();
G = basis.G(z);
Gz = basis.uz(z);
F = basis.F(z);
% All study numerical EVPs use the solved G formulation. Analytical bases
% expose the same normalized G/uz/F methods and do not require a solver.
if isempty(factors), error('IM21:InvalidNormalization','Missing normalization.'); end
end

function result = physicalResidual(basis,requestedKappa,c)
z = linspace(-c.D,0,c.referenceGridCount).';
w = [0.5;ones(numel(z)-2,1);0.5]/(numel(z)-1);
factors = basis.normalizationFactors();
G = basis.G(z);
Gz = basis.uz(z);
Gzz = basis.solver.evaluatePhysicalDerivative(basis.nativeModes,z,2)./factors;
lambda = basis.eigenvalues;
r = (basis.N2(z)-c.f0^2)/c.g;
terms = {-c.D^2*Gzz,(requestedKappa*c.D)^2*G,c.D^2*r.*G.*lambda};
residual = terms{1}+terms{2}-terms{3};
denominator = vecnorm(sqrt(w).*terms{1})+vecnorm(sqrt(w).*terms{2})+vecnorm(sqrt(w).*terms{3});
equation = vecnorm(sqrt(w).*residual)./max(denominator,realmin);
boundaries = [basis.evp.bottomBoundary,basis.evp.surfaceBoundary];
indices = [1,numel(z)]; boundaryErrors = zeros(2,numel(lambda));
scaleG = sqrt(sum(w.*(abs(G).^2+c.D^2*abs(Gz).^2),1));
for j = 1:2
    bc = boundaries(j); index = indices(j);
    boundaryResidual = -bc.a*G(index,:)+bc.b*Gz(index,:)-lambda.*(bc.c*G(index,:)-bc.d*Gz(index,:));
    dualScale = abs(bc.a)+abs(bc.b)/c.D+abs(lambda)*(abs(bc.c)+abs(bc.d)/c.D);
    boundaryErrors(j,:) = abs(boundaryResidual)./max(scaleG.*dualScale,realmin);
end
% This discretized-pencil backward residual is separate from off-grid
% physical-equation and boundary diagnostics; mixed row scaling is explicit.
[A,B] = requestedPencil(basis,requestedKappa);
V = basis.nativeModes;
backward = vecnorm(A*V-B*V.*lambda)./max((norm(A,2)+abs(lambda)*norm(B,2)).*vecnorm(V),realmin);
result = struct('equationError',max(equation(1:c.bands(end))),'boundaryError',max(boundaryErrors(:,1:c.bands(end)),[],'all'),'matrixBackwardError',max(backward(1:c.bands(end))), ...
    'equationPerMode',equation,'boundaryPerMode',boundaryErrors,'matrixBackwardPerMode',backward,'normalization','dimensionless D-scaled G equation and boundary dual trace scale; experimental');
% A concentration statistic, not an automatic localized-mode label.
energy = w.*(abs(G).^2+c.D^2*abs(Gz).^2);
result.upperQuarterFraction = sum(energy(z > -0.25*c.D,:),1)./sum(energy,1);
result.participationFraction = (sum(energy,1).^2)./(numel(z)*sum(energy.^2,1));
end

function [A,B] = requestedPencil(basis,kappa)
evp = IMInternalModes.waveModesAtWavenumber(N2=basis.N2,zDomain=basis.zDomain,f0=basis.evp.f0,g=basis.evp.g,k=kappa,surfaceBoundary=basis.evp.surfaceBoundary,bottomBoundary=basis.evp.bottomBoundary);
[A,B] = evp.assembleConfigured(basis.solver);
end

function ok = referencePasses(comparison,residual,c)
ok = all(comparison.eigenvalueError <= c.referenceEigenvalueTolerance) && all(comparison.frequencyError <= c.referenceEigenvalueTolerance) && all(comparison.subspaceError <= c.referenceSubspaceTolerance) && residual.equationError <= c.equationTolerance && residual.boundaryError <= c.boundaryTolerance;
end

function result = analyticalControl(problem,reference,requestedKappa,c)
result = struct('available',false,'reason','No exact analytical catalog for this profile.','comparison',[],'seconds',0);
if ~ismember(problem.profile,["constant","exponential"]), return; end
timer = tic;
try
    if problem.profile == "constant"
        solution = IMConstantStratificationSolution(N0=0.01,zDomain=problem.zDomain,f0=problem.f0,g=problem.g);
    else
        solution = IMExponentialStratificationSolution(N0=0.01,b=c.D,zDomain=problem.zDomain,f0=problem.f0,g=problem.g);
    end
    exact = solution.internalModes(reference.evp,nModes=c.referenceModes);
    result.comparison = compareBases(reference,exact,requestedKappa,c);
    result.available = true;
    result.reason = 'Exact analytical functions; normalization/integration and numerical root determination retain their own errors.';
catch exception
    result.reason = string(exception.identifier)+": "+string(exception.message);
end
result.seconds = toc(timer);
end

function result = failedComparison(c)
n = numel(c.bands);
result = struct('eigenvalueError',inf(1,n),'frequencyError',inf(1,n),'subspaceError',inf(1,n),'individualH1Error',inf(1,n),'FError',inf(1,n),'cutoffRelativeGap',nan(1,n),'minimumRelativeGap',NaN,'candidateLabels',[],'referenceLabels',[],'matching','unavailable');
end
function result = failedResidual()
result = struct('equationError',Inf,'boundaryError',Inf,'matrixBackwardError',Inf,'equationPerMode',[],'boundaryPerMode',[],'matrixBackwardPerMode',[],'normalization','unavailable','upperQuarterFraction',[],'participationFraction',[]);
end
function summary = finish(directory,records,configuration,seconds,complete)
summary = struct('completedCaseCount',numel(records),'requestedCaseCount',numel(configuration.selectedCaseIds),'elapsedSeconds',seconds,'caseMatrixComplete',complete,'productAssessmentStatus','deferred pending #19','productionQualification',false);
if ~isempty(records)
    summary.referenceQualifiedCount = nnz([records.referenceQualifiedExperimental]);
    summary.warmMatchesFullCount = nnz([records.warmMatchesFull]);
    summary.fallbackCount = nnz([records.fallbackRequired]);
    summary.medianWarmRawSeconds = median(arrayfun(@(r) r.warmCost.totalSeconds,records));
    summary.medianColdRawSeconds = median(arrayfun(@(r) r.coldCost.totalSeconds,records));
    summary.medianFullSolveSeconds = median([records.fullSolveSeconds]);
    summary.medianWarmWithValidationAndFallbackSeconds = median([records.warmTotalWithFallbackSeconds]);
    summary.medianWarmWithReferenceSeconds = median([records.validatedWithReferenceSeconds]);
    summary.medianNeighborAcquisitionSeconds = median([records.neighborAcquisitionSeconds]);
end
artifacts = dir(fullfile(directory,'*.mat'));
summary.retainedMatFileBytes = sum([artifacts.bytes]);
summary.retainedMatFiles = struct('name',{artifacts.name},'bytes',{artifacts.bytes});
save(fullfile(directory,'summary.mat'),'summary','records','configuration');
fid = fopen(fullfile(directory,'summary.json'),'w'); cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(summary,PrettyPrint=true));
disp(summary);
end
