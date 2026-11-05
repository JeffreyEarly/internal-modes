function report = validateNeutralLabels
% Check manuscript DD labels against independent analytical references.
%
% Deterministic controls use a shifted physical domain, separated bottom,
% interior and surface intervals, and a fully neutral column. Independent
% Gauss-Legendre quadrature and adaptive integrals check sine profiles,
% potential-enstrophy norms, complex nonmodal APV and displacement sources.
% Integration by parts includes both endpoint traces; affine sources cancel.
% The checks also cover signed and very small rotation, selected normalized
% common bases, repeated collection pages, empty sectors, and h=0 waves.
%
% Prescribed physical projection is checked separately from recovery. An
% intentionally unresolved quadrature must retain its physical error rather
% than acquiring a fitted sampled inverse. Error metrics are relative
% Frobenius errors for nonzero references and absolute Frobenius errors for
% zero references. Each comparison has an explicit tolerance in this file.
%
% Add this repository's tools folder to the MATLAB path and configure its
% declared package dependencies before calling. This tool temporarily adds
% its repository root and restores the caller's path and random-number state
% on completion or failure. Thread configuration remains caller controlled;
% the recorded audit used `maxNumCompThreads(1)` before calling this tool.
% No files are read or written and no research caches are required.
%
% ```matlab
% report = validateNeutralLabels;
% report.maximumErrors
% ```
%
% - Topic: Developer topics — Validate neutral labels
% - Declaration: report = validateNeutralLabels
% - Returns report: control count, maximum errors, runtime environment and repository provenance
% - Developer: true
arguments (Output)
    report (1,1) struct
end
callerPath = path;
callerRng = rng;
restoreState = onCleanup(@() restoreCallerState(callerPath,callerRng));
repo = fileparts(fileparts(mfilename("fullpath")));
addpath(repo);
domain = [-1337.5 25.5];
D = diff(domain);
intervals = [domain(1) -1198.25;-876.5 -421.25;-199.75 domain(2)];
counts = [3 4 2];
labels = repelem(1:3,counts);
numbers = [1:3 1:4 1:2];
lengths = diff(intervals,1,2);
m = pi*numbers./lengths(labels).';
amplitude = sqrt(2*D./lengths(labels)).';
[nodes,weights] = gaussRule(96);
z = [];
w = [];
for j = 1:3
    z = [z;intervals(j,1)+lengths(j)*(nodes+1)/2]; %#ok<AGROW>
    w = [w;lengths(j)*weights/2]; %#ok<AGROW>
end
metrics = struct();
checks = 0;
rotations = [7.2921e-5 7.2921e-5 1e-30 1e-100];
latitudes = [61 -37 42 -71];
for configuration = 1:numel(rotations)
    rotation = rotations(configuration);
    latitude = latitudes(configuration);
    f = 2*rotation*sind(latitude);
    basis = IMNeutralLabelModesBasis(intervals,zDomain=domain,nModesPerInterval=counts,rotationRate=rotation,latitude=latitude,g0=Inf,gd=-0.07);
    G = zeros(numel(z),sum(counts));
    Gz = G;
    for column = 1:sum(counts)
        rows = z > intervals(labels(column),1) & z < intervals(labels(column),2);
        phase = m(column)*(z(rows)-intervals(labels(column),1));
        G(rows,column) = amplitude(column)*sin(phase);
        Gz(rows,column) = amplitude(column)*m(column)*cos(phase);
    end
    Q = -f*Gz;
    verifyClose('sine',basis.G(z),G,2e-13);
    verifyClose('derivative',basis.Gz(z),Gz,2e-13);
    verifyClose('secondDerivative',basis.Gzz(z),-G.*m.^2,2e-13);
    verifyClose('APV',basis.APV(z),Q,2e-13);
    verifyClose('GNorm',G'*(w.*G)/D,eye(sum(counts)),2e-13);
    physicalGram = diag(f^2*D*m.^2);
    verifyClose('QNorm',Q'*(w.*Q),physicalGram,2e-13);
    verifyClose('reportedQNorm',basis.gramMatrix(variable="APV"),physicalGram,2e-13);
    endpointZ = sort(unique(intervals(:)));
    assert(all(basis.G(endpointZ)==0,'all')); checks = checks+1;
    assert(all(basis.Gzz(endpointZ)==0,'all')); checks = checks+1;
    projection = basis.projectionOnGrid(z,w,variable="APV",columns=[9 4 1 8]);
    coefficients = (1:9).'+1i*cos((1:9).');
    verifyClose('modalRecovery',projection.projectPhysical(Q*coefficients),coefficients([9 4 1 8]),1e-12);
    verifyClose('physicalDual',projection.projectPhysical(Q*coefficients),projection.projectPhysical(Q*coefficients),1e-15);
    assert(isempty(projection.recovery)); checks = checks+1;
    assert(projection.projectionKind=="prescribedDual"); checks = checks+1;
    assert(~projection.supportsGramAssessment && ~projection.supportsRoundTripAssessment); checks = checks+1;
    assert(isequal(projection.columnLabels,basis.columnLabels([9 4 1 8]))); checks = checks+1;
    verifyClose('selectedPhysicalGram',projection.provenance.continuousAPVGramMatrix,physicalGram([9 4 1 8],[9 4 1 8]),1e-15);
    % Arbitrary complex APV not lying in the retained cosine span.
    q = @(x) (0.7+0.2i)*exp((x-domain(1))/D)+sin(2.7*(x-domain(1))/D)+1i*((x-domain(1))/D).^4;
    sampled = f*q(z);
    expected = zeros(sum(counts),1);
    for column = 1:sum(counts)
        interval = intervals(labels(column),:);
        qz = @(x) amplitude(column)*m(column)*cos(m(column)*(x-interval(1)));
        expected(column) = -integral(@(x) qz(x).*q(x),interval(1),interval(2),AbsTol=1e-11,RelTol=1e-12)/(D*m(column)^2);
    end
    allProjection = basis.projectionOnGrid(z,w,variable="APV");
    verifyClose('nonmodalAdaptive',allProjection.projectPhysical(sampled),expected,4e-12);
    % Displacement source with nonzero traces: compare to the independently
    % integrated volume functional PLUS both endpoint traces.
    source = @(x) 1.2-0.3i+(0.7+0.1i)*(x-domain(1))/D+0.4*((x-domain(1))/D).^3+0.2i*sin((x-domain(1))/D);
    derivative = @(x) (0.7+0.1i)/D+1.2*((x-domain(1))/D).^2/D+0.2i*cos((x-domain(1))/D)/D;
    sourceExpected = zeros(sum(counts),1);
    volumeOnly = sourceExpected;
    for column = 1:sum(counts)
        interval = intervals(labels(column),:);
        profile = @(x) amplitude(column)*sin(m(column)*(x-interval(1)));
        volumeOnly(column) = integral(@(x) profile(x).*source(x),interval(1),interval(2),AbsTol=1e-11,RelTol=1e-12)/D;
        traces = amplitude(column)*m(column)*[-1 (-1)^numbers(column)].*source(interval);
        sourceExpected(column) = volumeOnly(column)+sum(traces)/(D*m(column)^2);
    end
    verifyClose('sourceIBP',allProjection.projectPhysical(-f*derivative(z)),sourceExpected,2e-12);
    assert(norm(volumeOnly-sourceExpected)>0.01); checks = checks+1;
    verifyClose('affineCancellation',allProjection.projectPhysical((-f*(0.3-0.2i)/D)*ones(size(z))),zeros(sum(counts),1),3e-13);
    % Generalized energy coefficients affect provenance only.
    changed = IMNeutralLabelModesBasis(intervals,zDomain=domain,nModesPerInterval=counts,rotationRate=rotation,latitude=latitude,g0=-Inf,gd=Inf);
    assert(isequal(changed.G(z),basis.G(z))); checks = checks+1;
    assert(isequal(changed.APV(z),basis.APV(z))); checks = checks+1;
    changedProjection = changed.projectionOnGrid(z,w,variable="APV");
    assert(isequal(changedProjection.physicalDualMatrix,allProjection.physicalDualMatrix)); checks = checks+1;
end
% Capable basis and repeated collection pages retain exact normalization,
% selected ordering, physical Gram and source projection semantics.
background = IMBackground(N2=@(x) neutralIntervalN2(x,intervals),zDomain=domain,latitude=-37,surfaceCondition="freeSurface");
raw = IMNeutralLabelModesBasis(intervals,zDomain=domain,nModesPerInterval=counts,latitude=-37);
rawProjection = raw.projectionOnGrid(z,w,variable="APV");
coefficients = [(1:9).'+2i,exp(1i*(1:9).')];
for normalization = ["rmsG","maxAbsG"]
    basisView = IMNeutralLabelModesBasis(background,neutralIntervals=intervals,nModesPerInterval=counts,normalization=normalization);
    if normalization == "rmsG", scale = ones(1,9); else, scale=amplitude; end
    verifyClose('basisViewProfiles',basisView.G(z),raw.G(z)./scale,1e-15);
    verifyClose('basisViewAPV',basisView.APV(z),raw.APV(z)./scale,1e-15);
    selected = basisView.selectColumns([9 2 5 1]);
    collection = IMBasisCollection({selected},basisIndex=[1 1 1]);
    subset = [3 1 4];
    for page = [3 1 3]
        actual = collection.projectionOnGrid(z,w,variable="APV",page=page,columns=subset);
        expectedMap = scale([5 9 1]).'.*rawProjection.physicalDualMatrix([5 9 1],:);
        verifyClose('collectionMap',actual.physicalDualMatrix,expectedMap,1e-15);
        expectedGram = raw.gramMatrix(variable="APV");
        expectedGram = expectedGram([5 9 1],[5 9 1])./(scale([5 9 1]).'*scale([5 9 1]));
        verifyClose('normalizedGram',actual.majorantGramMatrix,expectedGram,1e-15);
        verifyClose('provenanceGram',actual.provenance.continuousAPVGramMatrix,expectedGram,1e-15);
        verifyClose('basisViewCoefficients',actual.projectPhysical(raw.APV(z)*coefficients),scale([5 9 1]).'.*coefficients([5 9 1],:),1e-12);
        assert(isempty(actual.recovery)); checks=checks+1;
    end
    empty = basisView.selectColumns([]);
    emptyCollection = IMBasisCollection({empty},basisIndex=[1 1]);
    assert(isequal(size(emptyCollection.evaluate(z,variable="APV",pages=[2 1 2])),[numel(z) 0 3])); checks=checks+1;
    emptyProjection = emptyCollection.projectionOnGrid(z,w,variable="APV",page=2);
    assert(isequal(size(emptyProjection.physicalDualMatrix),[0 numel(z)])); checks=checks+1;
    assert(isequal(size(emptyProjection.projectPhysical(ones(numel(z),3))),[0 3])); checks=checks+1;
    assert(isequal(size(emptyProjection.projectPhysical(ones(numel(z),3))),[0 3])); checks=checks+1;
    assert(isempty(emptyProjection.provenance.continuousAPVGramMatrix)); checks=checks+1;
end
% Zero rotation retains geometry but correctly refuses APV inversion.
for opts = {struct('rotationRate',0),struct('latitude',0)}
    argumentsStruct = opts{1};
    argumentsStruct.zDomain = domain;
    argumentsStruct.nModesPerInterval = counts;
    pairs = namedargs2cell(argumentsStruct);
    zero = IMNeutralLabelModesBasis(intervals,pairs{:});
    assert(isequal(zero.APV(z),zeros(numel(z),9))); checks=checks+1;
    try
        zero.projectionOnGrid(z,w,variable="APV");
        error("IMNeutralLabelValidation:MissingZeroRotationFailure",'Expected zero-rotation failure.');
    catch exception
        assert(strcmp(exception.identifier,'IMNeutralLabelModesBasis:ZeroRotationProjection')); checks=checks+1;
    end
end
% Fully neutral column has only the positive sine ladder (no constant mode).
full = IMNeutralLabelModesBasis(domain,zDomain=domain,nModesPerInterval=5);
assert(isequal(full.localModeNumber,1:5)); checks=checks+1;
verifyClose('fullColumnWavenumbers',full.m_n,pi*(1:5)/D,1e-15);
assert(all(full.G(domain.')==0,'all')); checks=checks+1;
zFull = domain(1)+D*(nodes+1)/2;
wFull = D*weights/2;
fullProjection = full.projectionOnGrid(zFull,wFull,variable="APV");
fullExpected = zeros(5,1);
fullSourceExpected = zeros(5,1);
for column = 1:5
    wavenumber = pi*column/D;
    profile = @(x) sqrt(2)*sin(wavenumber*(x-domain(1)));
    profileDerivative = @(x) sqrt(2)*wavenumber*cos(wavenumber*(x-domain(1)));
    fullExpected(column) = -integral(@(x) profileDerivative(x).*q(x),domain(1),domain(2),AbsTol=1e-10,RelTol=1e-12)/(D*wavenumber^2);
    volume = integral(@(x) profile(x).*source(x),domain(1),domain(2),AbsTol=1e-10,RelTol=1e-12)/D;
    boundary = sqrt(2)*wavenumber*[-1 (-1)^column].*source(domain);
    fullSourceExpected(column) = volume+sum(boundary)/(D*wavenumber^2);
end
verifyClose('fullNonmodalAPV',fullProjection.projectPhysical(full.f0*q(zFull)),fullExpected,4e-12);
verifyClose('fullSourceIBP',fullProjection.projectPhysical(-full.f0*derivative(zFull)),fullSourceExpected,2e-12);
verifyClose('fullAffineCancellation',fullProjection.projectPhysical((-full.f0*(0.3-0.2i)/D)*ones(size(zFull))),zeros(5,1),3e-13);
% A deliberately unresolved fixed quadrature must retain its physical error;
% the source dual may not turn into a fitted sample inverse.
coarseZ = domain(1)+D*[0.03;0.19;0.53;0.9];
coarseW = D*[0.2;0.2;0.2;0.4];
coarseQ = -full.f0*sqrt(2)*(pi*(1:5)/D).*cos((coarseZ-domain(1))*(pi*(1:5)/D));
coarse = full.projectionOnGrid(coarseZ,coarseW,variable="APV");
coarseReference = diag(1./(full.f0^2*D*(pi*(1:5)/D).^2))*(coarseQ'.*coarseW.');
verifyClose('unresolvedPhysicalMap',coarse.physicalDualMatrix,coarseReference,4e-15);
assert(norm(coarse.physicalDualMatrix*coarseQ-eye(5),'fro')>0.1); checks=checks+1;
verifyClose('unresolvedPhysicalApplication',coarse.projectPhysical(coarseQ),coarse.projectPhysical(coarseQ),1e-15);
% The same DD auxiliary shape underlies zero-depth waves. Check their
% independent wave-energy normalization rather than confusing labels/h=0.
critical = IMBackground(N2=@(x) background.f0^2*ones(size(x)),zDomain=domain,latitude=-37,surfaceCondition="freeSurface");
k = 3.7e-4;
wave = IMZeroDepthWaveModesBasis(critical,k,criticalIntervals=domain,nModesPerInterval=5);
mFull = pi*(1:5)/D;
phase = (zFull-domain(1))*mFull;
expectedG = sqrt(2)*sin(phase)./(D*sqrt(mFull.^2+k^2));
expectedF = sqrt(2)*mFull.*cos(phase)./sqrt(mFull.^2+k^2);
verifyClose('zeroDepthG',wave.G(zFull),expectedG,3e-13);
verifyClose('zeroDepthF',wave.F(zFull),expectedF,3e-13);
verifyClose('zeroDepthWaveEnergy',sum(wFull.*(abs(wave.F(zFull)).^2+(k*D)^2*abs(wave.G(zFull)).^2),1)/D,ones(1,5),3e-13);
assert(isequal(wave.h,zeros(1,5))); checks=checks+1;
assert(isequal(wave.frequency,abs(critical.f0)*ones(1,5))); checks=checks+1;
report = struct(checks=checks,maximumErrors=metrics,matlabVersion=string(version),computer=string(computer),repositoryRoot=string(repo));

    function verifyClose(name,actual,expected,tolerance)
        assert(isequal(size(actual),size(expected)),"IMNeutralLabelValidation:ShapeMismatch",'%s shape mismatch.',name);
        if norm(expected,'fro')==0
            errorValue=norm(actual,'fro');
        else
            errorValue=norm(actual-expected,'fro')/norm(expected,'fro');
        end
        if ~isfield(metrics,name), metrics.(name)=0; end
        metrics.(name)=max(metrics.(name),errorValue);
        assert(isfinite(errorValue) && errorValue<tolerance,"IMNeutralLabelValidation:ReferenceMismatch",'%s relative error %.17g exceeds %.17g.',name,errorValue,tolerance);
        checks=checks+1;
    end
end

function [nodes,weights] = gaussRule(n)
index=(1:n-1).';
offDiagonal=index./sqrt(4*index.^2-1);
[vectors,values]=eig(diag(offDiagonal,1)+diag(offDiagonal,-1));
[nodes,order]=sort(diag(values));
weights=2*vectors(1,order).'.^2;
end

function restoreCallerState(callerPath,callerRng)
path(callerPath);
rng(callerRng);
end

function N2 = neutralIntervalN2(z,intervals)
N2 = 1e-4*ones(size(z));
for interval = 1:size(intervals,1)
    N2(z >= intervals(interval,1) & z <= intervals(interval,2)) = 0;
end
end
