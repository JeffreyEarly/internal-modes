function report = runNeutralGeostrophicValidation(options)
% Qualify balanced neutral-cap references before a coordinate-study rerun.
arguments
    options.nModes (1,1) double {mustBeInteger,mustBePositive} = 40
    options.resolutions (1,2) double {mustBeInteger,mustBePositive} = [384 512]
    options.retryResolution (1,1) double {mustBeInteger,mustBePositive} = 1024
    options.discretization (1,1) string {mustBeMember(options.discretization,["global","piecewise"])} = "global"
    options.caseNames (1,:) string = strings(1,0)
    options.outputPath (1,1) string = ""
end
repoRoot = fileparts(fileparts(fileparts(mfilename("fullpath"))));
oldPath = path;
cleanup = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(repoRoot,fullfile(repoRoot,"UnitTestsV2","TestSupport"));
if options.outputPath == ""
    options.outputPath = fullfile(fileparts(mfilename("fullpath")),"readiness-report.json");
end
profile = neutralGeostrophicProfile();
cases = validationCases(profile);
availableNames = strings(1,numel(cases));
for i = 1:numel(cases)
    availableNames(i) = cases{i}.name;
end
if ~isempty(options.caseNames)
    unknown = setdiff(options.caseNames,availableNames);
    if ~isempty(unknown)
        error("NeutralGeostrophicValidation:UnknownCase","Unknown validation case(s): %s.",join(unknown,", "));
    end
    selected = ismember(availableNames,options.caseNames);
    cases = cases(selected);
    availableNames = availableNames(selected);
end
report.schemaVersion = 2;
report.startedAt = string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ssXXX"));
report.profile = rmfield(profile,"N2");
report.settings = struct("nIndexedModes",options.nModes,"includeRequiredSignedAndConstantSectors",true, ...
    "coordinateKind","z","resolutions",options.resolutions,"retryResolution",options.retryResolution, ...
    "shootingTolerances",[1e-11 1e-13],"quadraturePointsPerInterval",[257 513 1025], ...
    "referenceTolerance",1e-7,"quadratureTolerance",1e-9,"gramTolerance",2e-7, ...
    "rotationRate",7.2921e-5,"latitude",45,"g",9.81,"comparisonStudyRerun",false, ...
    "caseNames",availableNames);
report.settings.discretization = options.discretization;
report.settings.interfaceJoins = profile.joins(2:end-1);
report.settings.resolutionMeans = "Total coefficients per solved field, summed across intervals; diagnostic antiderivatives may have higher degree.";
report.sourcesAtStart = sourceIdentities(repoRoot);
report.cases = cell(1,numel(cases));
report.readyForNeutralBalancedComparisons = false;
report.status = "running";
for i = 1:numel(cases)
    configuration = cases{i};
    fprintf("Neutral readiness %d/%d: %s\n",i,numel(cases),configuration.name);
    timer = tic;
    row = struct("name",configuration.name,"family",configuration.family,"parameters",configuration.parameters, ...
        "qualification","inconclusive","reason","","seconds",NaN);
    try
        problem = configuration.create();
        coarse = neutralGeostrophicReference(problem,nModes=options.nModes,includeSignedSectors=true,tolerance=1e-11);
        tight = neutralGeostrophicReference(problem,nModes=numel(coarse.h),rootBrackets=coarse.rootBrackets,tolerance=1e-13);
        row.referenceScan = coarse.scan;
        row.modeNumber = tight.modeNumber;
        row.inverseDepth = tight.inverseDepth;
        row.infiniteEquivalentDepth = isinf(tight.h);
        row.equivalentDepth = tight.h;
        row.publicEigenvalues = tight.eigenvalues;
        [z,weights] = neutralGeostrophicQuadrature(profile.joins,513);
        row.independentReferenceAgreement = compactMetrics(neutralGeostrophicComparison(tight,coarse,z,weights));
        requested = numel(tight.h);
        [first,firstFailure] = attemptSolve(problem,options.resolutions(1),requested,options.discretization,profile.joins);
        [second,secondFailure] = attemptSolve(problem,options.resolutions(2),requested,options.discretization,profile.joins);
        row.spectralSolveFailures = [firstFailure secondFailure];
        if ~isempty(first) && ~isempty(second)
            row.spectralRefinement = compactMetrics(compareBases(tight,first,second,z,weights));
        else
            row.spectralRefinement = struct("maximumError",Inf,"reason","One or both initial spectral references failed qualification.");
        end
        row.finalResolution = options.resolutions(2);
        candidate = second;
        if row.spectralRefinement.maximumError > 1e-7
            row.finalResolution = options.retryResolution;
            [candidate,retryFailure] = attemptSolve(problem,options.retryResolution,requested,options.discretization,profile.joins);
            row.spectralSolveFailures = [row.spectralSolveFailures retryFailure];
            if isempty(candidate)
                error("NeutralGeostrophicValidation:SpectralReferenceFailure", ...
                    "The spectral reference failed qualification at the retry resolution %d.",options.retryResolution);
            end
            if ~isempty(second)
                row.spectralRefinement = compactMetrics(compareBases(tight,second,candidate,z,weights));
            end
        end
        row.dominantGOrientationFallback = false(1,requested);
        if isfield(candidate.metadata,"nPerInterval")
            row.pointsPerInterval = candidate.metadata.nPerInterval;
        end
        if isfield(candidate.metadata,"coefficientCleanup")
            row.coefficientCleanup = candidate.metadata.coefficientCleanup;
        end
        if isfield(candidate.metadata,"dominantGOrientationFallback")
            row.dominantGOrientationFallback = candidate.metadata.dominantGOrientationFallback;
        end
        [z1,w1] = neutralGeostrophicQuadrature(profile.joins,257);
        metrics1 = neutralGeostrophicComparison(tight,candidate,z1,w1);
        metrics2 = neutralGeostrophicComparison(tight,candidate,z,weights);
        row.quadratureConvergence = metricDistance(metrics1,metrics2);
        row.finalQuadraturePointsPerInterval = 513;
        if row.quadratureConvergence > 1e-9
            [z,weights] = neutralGeostrophicQuadrature(profile.joins,1025);
            metrics3 = neutralGeostrophicComparison(tight,candidate,z,weights);
            row.quadratureConvergence = metricDistance(metrics2,metrics3);
            row.finalQuadraturePointsPerInterval = 1025;
            metrics2 = metrics3;
        end
        row.productionAgainstIndependent = compactMetrics(metrics2);
        row.mathematicalChecks = mathematicalChecks(tight,candidate,z,weights);
        checks = [row.independentReferenceAgreement.maximumError <= 1e-7, ...
            row.spectralRefinement.maximumError <= 1e-7,metrics2.maximumError <= 1e-7, ...
            row.quadratureConvergence <= 1e-9,row.mathematicalChecks.maximumResidual <= 1e-7, ...
            row.mathematicalChecks.maximumGramError <= 2e-7];
        row.acceptanceChecks = struct("independentReference",checks(1),"spectralRefinement",checks(2), ...
            "productionAccuracy",checks(3),"quadratureConvergence",checks(4), ...
            "equationsAndBoundaries",checks(5),"gramIdentities",checks(6));
        if all(checks)
            row.qualification = "qualified";
        else
            row.reason = "At least one required accuracy or reference-convergence gate remains unresolved.";
        end
    catch exception
        row.reason = string(exception.message);
        row.errorIdentifier = string(exception.identifier);
        row.errorStack = string({exception.stack.name});
    end
    row.seconds = toc(timer);
    report.cases{i} = row;
    writeReport(options.outputPath,report);
    fprintf("  %s (%.1f s)\n",row.qualification,row.seconds);
end
qualified = cellfun(@(row) row.qualification == "qualified",report.cases);
report.readyForNeutralBalancedComparisons = all(qualified);
report.qualifiedCaseCount = nnz(qualified);
report.caseCount = numel(cases);
report.status = "complete";
report.finishedAt = string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ssXXX"));
report.sourcesAtFinish = sourceIdentities(repoRoot);
report.sourcesChangedDuringValidation = ~isequal(report.sourcesAtStart,report.sourcesAtFinish);
if report.sourcesChangedDuringValidation
    report.readyForNeutralBalancedComparisons = false;
    report.sourceQualification = "Source changes during validation require review or rerunning the affected cases.";
else
    report.sourceQualification = "Source identities remained fixed during validation.";
end
writeReport(options.outputPath,report);
end

function cases = validationCases(p)
rotation = 7.2921e-5;
f0 = 2*rotation*sind(45);
cases = cell(1,0);
boundaries = {"manuscript",0,-0.01;"wvm",-p.integratedN2,p.integratedN2};
for i = 1:size(boundaries,1)
    name = boundaries{i,1};
    g0 = boundaries{i,2};
    gd = boundaries{i,3};
    parameters = struct("g0",g0,"gd",gd,"surfaceBoundary","freeSurface");
    cases{end+1} = struct("name","apv-"+name,"family","apv","parameters",parameters, ...
        "create",@() IMInternalModes.geostrophicAPVModes(N2=p.N2,zDomain=p.zDomain,rotationRate=rotation,latitude=45,g0=g0,gd=gd)); %#ok<AGROW>
    cases{end+1} = struct("name","mda-"+name,"family","mda","parameters",parameters, ...
        "create",@() IMInternalModes.meanDensityAnomalyModes(N2=p.N2,zDomain=p.zDomain,rotationRate=rotation,latitude=45,g0=g0,gd=gd)); %#ok<AGROW>
end
for wavelength = [10e3 1e3]
    k = 2*pi/wavelength;
    parameters = struct("wavelength",wavelength,"k",k,"alpha0","stratification default","alphaD","Inf");
    cases{end+1} = struct("name","gpe-default-"+wavelength,"family","gpe","parameters",parameters, ...
        "create",@() IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(N2=p.N2,zDomain=p.zDomain,k=k,rotationRate=rotation,latitude=45)); %#ok<AGROW>
    alpha0 = 2*f0^2/diff(p.zDomain);
    alphaD = 3*f0^2/diff(p.zDomain);
    parameters.alpha0 = alpha0;
    parameters.alphaD = alphaD;
    cases{end+1} = struct("name","gpe-finite-"+wavelength,"family","gpe","parameters",parameters, ...
        "create",@() IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(N2=p.N2,zDomain=p.zDomain,k=k,rotationRate=rotation,latitude=45,alpha0=alpha0,alphaD=alphaD)); %#ok<AGROW>
    parameters.alpha0 = "Inf";
    parameters.alphaD = "Inf";
    cases{end+1} = struct("name","gpe-inactive-"+wavelength,"family","gpe","parameters",parameters, ...
        "create",@() IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(N2=p.N2,zDomain=p.zDomain,k=k,rotationRate=rotation,latitude=45,alpha0=Inf,alphaD=Inf)); %#ok<AGROW>
    parameters = struct("wavelength",wavelength,"k",k,"surfaceBoundary","freeSurface","endpointResponses",["surface","bottom"]);
    cases{end+1} = struct("name","zeroAPV-"+wavelength,"family","zeroAPV","parameters",parameters, ...
        "create",@() IMGeostrophicZeroAPVModes(N2=p.N2,zDomain=p.zDomain,k=k,rotationRate=rotation,latitude=45)); %#ok<AGROW>
end
end

function basis = solveProblem(problem,nEVP,nModes,discretization,joins)
if discretization == "piecewise"
    basis = neutralGeostrophicPiecewiseBasis(problem,nEVP,nModes=nModes,joins=joins);
    return
end
solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="z");
if isa(problem,"IMGeostrophicZeroAPVModes")
    basis = solver.solveGeostrophicZeroAPVModes(problem);
else
    basis = solver.solveEVP(problem,nModes=nModes);
end
end

function [basis,failure] = attemptSolve(problem,nEVP,nModes,discretization,joins)
failure = struct("resolution",{},"identifier",{},"message",{});
try
    basis = solveProblem(problem,nEVP,nModes,discretization,joins);
catch exception
    basis = [];
    failure = struct("resolution",nEVP,"identifier",string(exception.identifier),"message",string(exception.message));
end
end

function metrics = compareBases(reference,first,second,z,weights)
firstReference = reference;
firstReference.evaluate = @(points) sampleNeutralGeostrophicBasis(first,points);
if ~isa(first,"IMGeostrophicZeroAPVModesBasis")
    firstReference.h = first.h;
    firstReference.inverseDepth = 1./first.h;
    firstReference.modeNumber = first.modeNumber;
end
metrics = neutralGeostrophicComparison(firstReference,second,z,weights);
end

function distance = metricDistance(first,second)
if first.matchedModeCount ~= second.matchedModeCount || first.matchedModeCount == 0
    distance = Inf;
    return
end
distance = max(abs([first.FH1Error-second.FH1Error first.GH1Error-second.GH1Error]));
for name = ["sampledFGram","sampledGGram"]
    delta = norm(first.(name)-second.(name),"fro")/max(norm(second.(name),"fro"),realmin);
    distance = max(distance,delta);
end
end

function metrics = compactMetrics(metrics)
metrics = rmfield(metrics,["sampledFGram","sampledGGram"]);
end

function result = mathematicalChecks(reference,basis,z,weights)
c = reference.configuration;
v = sampleNeutralGeostrophicBasis(basis,z);
D = diff(c.zDomain);
r = reference.inverseDepth;
FScale = max(abs(v.F),[],1);
GScale = max(abs(v.G),[],1);
first = v.Fz+(c.N2(z)/c.g).*v.G;
firstScale = max(FScale/D,max(abs((c.N2(z)/c.g).*v.G),[],1));
firstError = max(abs(first),[],1)./max(firstScale,realmin);
if c.family == "mda"
    % MDA's diagnostic pressure differs by its surface constant.
    if isstruct(basis)
        derivative = v.Gzz;
    else
        derivative = basis.solver.evaluatePhysicalDerivative(basis.nativeModes,z,2)./basis.normalizationFactors(basis.normalization);
    end
    second = derivative+r.*(c.N2(z)/c.g).*v.G;
    secondScale = max(GScale/D^2,max(abs(r.*(c.N2(z)/c.g).*v.G),[],1));
else
    second = v.Gz-r.*v.F;
    secondScale = max(GScale/D,max(abs(r.*v.F),[],1));
end
secondError = max(abs(second),[],1)./max(secondScale,realmin);
FGram = v.F.'*(weights.*v.F);
GGram = v.G.'*(weights.*(c.N2(z)/c.g).*v.G);
Fs = v.F(end,:);
Fb = v.F(1,:);
Gs = v.G(end,:);
Gb = v.G(1,:);
if c.family == "apv"
    if c.g0 ~= 0
        GGram = GGram+(c.sigma+c.g/c.g0)*(Fs.'*Fs);
    end
    if c.gd ~= 0
        GGram = GGram+(c.g/c.gd)*(Fb.'*Fb);
    end
    FGramError = norm(FGram/D-eye(numel(r)),2);
    scale = sqrt(max(abs(r),1/D));
    GGramError = norm((GGram/D-diag(r))./(scale.'*scale),2);
elseif c.family == "gpe"
    FGram = FGram+(c.f0^2/c.alpha0)*(Fs.'*Fs)+(c.f0^2/c.alphaD)*(Fb.'*Fb);
    GGram = GGram+(1-c.g*c.k^2/c.alpha0)*(Fs.'*Fs)-(c.g*c.k^2/c.alphaD)*(Fb.'*Fb);
    FGramError = norm(FGram/D-eye(numel(r)),2);
    scale = sqrt(max(abs(r),1/D));
    GGramError = norm((GGram/D-diag(r))./(scale.'*scale),2);
elseif c.family == "mda"
    if isfinite(c.g0)
        GGram = GGram+(c.g0/c.g)*(Gs.'*Gs);
    end
    if isfinite(c.gd)
        GGram = GGram+(c.gd/c.g)*(Gb.'*Gb);
    end
    FGramError = 0;
    GGramError = norm(GGram-diag(basis.signatures),2);
else
    endpoint = [Gs-c.sigma*Fs;Gb];
    FGramError = norm(endpoint-eye(2),2);
    GGramError = 0;
end
boundaryError = endpointErrors(c,v,r,reference.eigenvalues,FScale,GScale,D);
result = struct("firstEquationError",firstError,"secondEquationError",secondError, ...
    "boundaryError",boundaryError,"maximumResidual",max([firstError secondError boundaryError]),"FGramError",FGramError, ...
    "GGramError",GGramError,"maximumGramError",max(FGramError,GGramError));
result.interfaceError = 0;
if isfield(basis.metadata,"interfaceResiduals")
    result.interfaceError = max(basis.metadata.interfaceResiduals,[],"all");
    result.maximumResidual = max(result.maximumResidual,result.interfaceError);
end
end

function errors = endpointErrors(c,v,r,Lambda,FScale,GScale,D)
Fs = v.F(end,:);
Fb = v.F(1,:);
Gs = v.G(end,:);
Gb = v.G(1,:);
if c.family == "apv"
    if isinf(c.g0)
        surface = Gs-c.sigma*Fs;
        surfaceScale = GScale+c.sigma*FScale;
    else
        surface = c.g0*Gs-(c.g+c.sigma*c.g0)*Fs;
        surfaceScale = abs(c.g0)*GScale+abs(c.g+c.sigma*c.g0)*FScale;
    end
    if isinf(c.gd)
        bottom = Gb;
        bottomScale = GScale;
    else
        bottom = c.g*Fb+c.gd*Gb;
        bottomScale = c.g*FScale+abs(c.gd)*GScale;
    end
elseif c.family == "gpe"
    if isinf(c.alpha0)
        surface = Gs-Fs;
        surfaceScale = GScale+FScale;
    else
        surface = c.alpha0*(Gs-Fs)+c.g*Lambda.*Fs;
        surfaceScale = c.alpha0*(GScale+FScale)+c.g*abs(Lambda).*FScale;
    end
    if isinf(c.alphaD)
        bottom = Gb;
        bottomScale = GScale;
    else
        bottom = c.alphaD*Gb-c.g*Lambda.*Fb;
        bottomScale = c.alphaD*GScale+c.g*abs(Lambda).*FScale;
    end
elseif c.family == "mda"
    if isinf(c.g0)
        surface = Gs;
        surfaceScale = GScale;
    else
        surface = v.Gz(end,:)-(c.g0/c.g)*r.*Gs;
        surfaceScale = GScale/D+abs(c.g0*r/c.g).*GScale;
    end
    if isinf(c.gd)
        bottom = Gb;
        bottomScale = GScale;
    else
        bottom = v.Gz(1,:)+(c.gd/c.g)*r.*Gb;
        bottomScale = GScale/D+abs(c.gd*r/c.g).*GScale;
    end
else
    surface = Gs-c.sigma*Fs-[1 0];
    bottom = Gb-[0 1];
    surfaceScale = max(1,GScale+c.sigma*FScale);
    bottomScale = max(1,GScale);
end
errors = [abs(surface)./max(surfaceScale,realmin) abs(bottom)./max(bottomScale,realmin)];
end

function identities = sourceIdentities(repoRoot)
names = ["@IMGeostrophicModes","@IMGeostrophicModesBasis","@IMEigenvalueProblem","@IMBasisSet", ...
    "@IMSolver","@IMSolverSpectral","@IMSolverFiniteDifference","@IMGeostrophicFiniteDifferenceTools", ...
    "@IMInternalModes","@IMInternalModesBasis", ...
    "@IMGeostrophicZeroAPVModes","@IMGeostrophicZeroAPVModesBasis", ...
    "@IMMeanDensityAnomalyModes","@IMMeanDensityAnomalyModesBasis","@IMGeostrophicTransform", ...
    "UnitTestsV2/TestSupport","tools/neutral-geostrophic-validation"];
identities = struct("path",{},"sha256",{});
for folder = names
    files = dir(fullfile(repoRoot,folder,"*.m"));
    for i = 1:numel(files)
        file = fullfile(files(i).folder,files(i).name);
        bytes = unicode2native(fileread(file),"UTF-8");
        hash = java.security.MessageDigest.getInstance("SHA-256");
        hash.update(bytes);
        digest = typecast(hash.digest(),"uint8");
        identities(end+1) = struct("path",string(erase(file,[repoRoot filesep])), ...
            "sha256",string(lower(reshape(dec2hex(digest,2).',1,[])))); %#ok<AGROW>
    end
end
files = ["@InternalModesSpectral/InternalModesSpectral.m", ...
    "IMGeostrophicAPVEndpointInertia.m","IMInternalModesNormalizationTools.m","IMModeOrientationTools.m","Normalization.m", ...
    "IMConstantGeostrophicAPVCatalog.m","IMExponentialGeostrophicAPVCatalog.m"];
for file = files
    fullPath = fullfile(repoRoot,file);
    if isfile(fullPath)
        hash = java.security.MessageDigest.getInstance("SHA-256");
        hash.update(unicode2native(fileread(fullPath),"UTF-8"));
        digest = typecast(hash.digest(),"uint8");
        identities(end+1) = struct("path",file,"sha256",string(lower(reshape(dec2hex(digest,2).',1,[])))); %#ok<AGROW>
    end
end
dependencyPath = which("standardChop");
if ~isempty(dependencyPath)
    hash = java.security.MessageDigest.getInstance("SHA-256");
    hash.update(unicode2native(fileread(dependencyPath),"UTF-8"));
    digest = typecast(hash.digest(),"uint8");
    identities(end+1) = struct("path",string(dependencyPath),"sha256",string(lower(reshape(dec2hex(digest,2).',1,[]))));
end
end

function writeReport(file,report)
fid = fopen(file,"w");
if fid < 0
    error("NeutralGeostrophicValidation:CannotWriteReport","Cannot write %s.",file);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,"%s\n",jsonencode(report,PrettyPrint=true));
end
