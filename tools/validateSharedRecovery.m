function report = validateSharedRecovery
% Check shared recovery against independent complex numerical references.
%
% The deterministic controls compare least squares with a reference SVD,
% corrected physical recovery with $$(BA)^{-1}B$$, and noise gain with the
% smallest singular value of $$DAC$$. They also check observation units,
% extended-span selection, changes of physical test basis and extreme scales.
% The reference pseudoinverse is used only here, on full-rank test matrices.
% No additional package dependencies are required.
%
% The caller's MATLAB path and random-number state are restored on completion
% or failure. Thread configuration remains under caller control; for the
% recorded qualification run, call `maxNumCompThreads(1)` before this tool.
%
% - Topic: Developer topics — Validate recovery
% - Declaration: report = validateSharedRecovery
% - Returns report: trial counts, maximum relative errors and extreme-scale results
% - Developer: true
arguments (Output)
    report (1,1) struct
end
callerPath = path;
callerRng = rng;
restoreState = onCleanup(@() restoreCallerState(callerPath,callerRng));
repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repositoryRoot);
rng(82731,"twister");
tolerance = 1e-10;
errors = struct(leastSquares=0,correctedPhysicalDual=0,observationUnits=0,observationNoiseGain=0,extendedSpan=0);
for trial = 1:40
    m = randi([4 10]);
    n = randi([1 m-1]);
    A = randn(m,n)+1i*randn(m,n);
    B = randn(n,m)+1i*randn(n,m);
    C = randn(n,n)+1i*randn(n,n)+2*eye(n);
    D = randn(m,m)+1i*randn(m,m)+4*eye(m);
    F = randn(m+2,m)+1i*randn(m+2,m);
    fitted = prepareIMRecovery(A,coefficientScale=C,observationFactor=D,fitFactor=F);
    expected = pinv(F*A)*F;
    errors.leastSquares = max(errors.leastSquares,relativeError(fitted.forwardMatrix,expected));
    physical = prepareIMRecovery(A,policy="correctedPhysicalDual",coefficientScale=C,observationFactor=D,physicalDual=B);
    expected = (B*A)\B;
    errors.correctedPhysicalDual = max(errors.correctedPhysicalDual,relativeError(physical.forwardMatrix,expected));
    U = diag(10.^linspace(-5,5,m).*exp(1i*randn(1,m)));
    units = prepareIMRecovery(U*A,policy="correctedPhysicalDual",coefficientScale=C,observationFactor=D/U,physicalDual=B/U);
    errors.observationUnits = max(errors.observationUnits,relativeError(units.forwardMatrix*U,physical.forwardMatrix));
    optimal = prepareIMRecovery(A,coefficientScale=C,observationFactor=D);
    noiseGain = norm((C\optimal.forwardMatrix)/D,2);
    expectedNoiseGain = 1/min(svd(D*A*C));
    errors.observationNoiseGain = max(errors.observationNoiseGain,abs(noiseGain-expectedNoiseGain)/noiseGain);
    extended = prepareIMRecovery(A(:,1),policy="extendedSpan",columnLabels="one",observationFactor=D,extendedSynthesis=A,extendedColumnLabels=["one" string(2:n)],extendedActiveColumnMask=true(1,n),retainedColumnIndices=1,extendedCoefficientScale=C);
    errors.extendedSpan = max(errors.extendedSpan,relativeError(extended.forwardMatrix,optimal.forwardMatrix(1,:)));
end

rng(927,"twister");
m = 7;
n = 3;
A = randn(m,n)+1i*randn(m,n);
B = randn(n,m)+1i*randn(n,m);
physical = prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=B,observationFactor=eye(m));
errors.physicalTestBasis = 0;
for trial = 1:30
    L = randn(n,n)+1i*randn(n,n)+2*eye(n);
    transformed = prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=L*B,observationFactor=eye(m));
    errors.physicalTestBasis = max(errors.physicalTestBasis,relativeError(transformed.forwardMatrix,physical.forwardMatrix));
end

names = ["huge dual";"huge synthesis and dual";"extreme physical row units";"large physical row norm";"subnormal raw overlap";"tiny synthesis";"tiny dual";"overflowing physical row norm";"overflowing synthesis norm"];
factories = {@() prepareIMRecovery(eye(2),policy="correctedPhysicalDual",physicalDual=1e200*eye(2),observationFactor=eye(2)); ...
    @() prepareIMRecovery(1e200*eye(2),policy="correctedPhysicalDual",physicalDual=1e200*eye(2),observationFactor=eye(2)); ...
    @() prepareIMRecovery(eye(2),policy="correctedPhysicalDual",physicalDual=diag([1e300 1e-300]),observationFactor=eye(2)); ...
    @() prepareIMRecovery([1;1],policy="correctedPhysicalDual",physicalDual=[1e308 1e308],observationFactor=eye(2)); ...
    @() prepareIMRecovery(1e-200*eye(2),policy="correctedPhysicalDual",physicalDual=1e-200*eye(2),observationFactor=eye(2)); ...
    @() prepareIMRecovery(1e-200*eye(2),observationFactor=eye(2)); ...
    @() prepareIMRecovery(eye(2),policy="correctedPhysicalDual",physicalDual=1e-300*eye(2),observationFactor=eye(2)); ...
    @() prepareIMRecovery([1;0;0;0],policy="correctedPhysicalDual",physicalDual=1e308*ones(1,4),observationFactor=eye(4)); ...
    @() prepareIMRecovery(1e308*ones(4,1),observationFactor=eye(4))};
expectedMaps = {eye(2);1e-200*eye(2);eye(2);0.5*ones(1,2);1e200*eye(2);1e200*eye(2);eye(2);ones(1,4);2.5e-309*ones(1,4)};
relativeErrors = zeros(numel(names),1);
roundTripErrors = zeros(numel(names),1);
for index = 1:numel(names)
    prepared = factories{index}();
    assert(prepared.success,"IMRecoveryValidation:UnexpectedFailure","Extreme-scale control '%s' failed: %s.",names(index),prepared.status);
    relativeErrors(index) = relativeError(prepared.forwardMatrix,expectedMaps{index});
    roundTripErrors(index) = prepared.diagnostics.roundTripError;
end
errors.extremeScales = max(relativeErrors);

errors.scaledWeakEquation = 0;
errors.scaledWeakMap = 0;
for alpha = [1 1e100 1e-100 1e200 1e-200]
    A = alpha*[1 0;0 1;1 1];
    B = alpha*[1 0 0;0 1 0];
    prepared = prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=B,observationFactor=eye(3));
    errors.scaledWeakMap = max(errors.scaledWeakMap,relativeError(alpha*prepared.forwardMatrix,[1 0 0;0 1 0]));
    errors.scaledWeakEquation = max(errors.scaledWeakEquation,prepared.diagnostics.scaledDefiningEquationResidual);
end

emptyObservationStatuses = strings(3,1);
counts = [0 1 3];
for index = 1:numel(counts)
    prepared = prepareIMRecovery(zeros(0,counts(index)),fitFactor=zeros(1,0),failureMode="diagnostic");
    emptyObservationStatuses(index) = prepared.status;
    if counts(index) == 0
        assert(prepared.success,"IMRecoveryValidation:EmptySpanFailure","The empty span should recover without observations.");
    else
        assert(~prepared.success && isempty(prepared.forwardMatrix),"IMRecoveryValidation:UnexpectedRecovery","An active span cannot recover without observations.");
    end
end
allErrors = cell2mat(struct2cell(errors));
assert(all(isfinite(allErrors)) && all(allErrors < tolerance),"IMRecoveryValidation:ReferenceMismatch","A shared recovery control exceeded relative tolerance %.3g; maximum error %.3g.",tolerance,max(allErrors));
report = struct(complexTrials=40,physicalTestBasisTrials=30,randomSeeds=[82731 927],relativeTolerance=tolerance,maxRelativeErrors=errors,extremeScales=table(names,relativeErrors,roundTripErrors),emptyObservationStatuses=emptyObservationStatuses,matlabVersion=string(version));
fprintf("Shared recovery audit passed: 40 complex controls, 30 physical test-basis controls, %d extreme scales.\n",numel(names));
fprintf("Maximum relative errors: LS %.3g; physical %.3g; units %.3g; noise %.3g; extended %.3g; test basis %.3g; extremes %.3g.\n",errors.leastSquares,errors.correctedPhysicalDual,errors.observationUnits,errors.observationNoiseGain,errors.extendedSpan,errors.physicalTestBasis,errors.extremeScales);
fprintf("Scaled weak controls: map error %.3g; defining-equation residual %.3g. Empty-observation controls passed.\n",errors.scaledWeakMap,errors.scaledWeakEquation);
end

function value = relativeError(actual,expected)
scale = max([abs(real(expected(:)));abs(imag(expected(:)));realmin]);
value = norm(actual/scale-expected/scale,2)/max(norm(expected/scale,2),realmin);
end

function restoreCallerState(callerPath,callerRng)
path(callerPath);
rng(callerRng);
end
