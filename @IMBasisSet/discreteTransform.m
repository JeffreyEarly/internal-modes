function [transform, assessment] = discreteTransform(self, options)
% Build and assess a scalar Galerkin transform on a fixed quadrature rule.
%
% A transform can be built from explicit physical points `z`, or from an
% exact requested point count `nPoints`. The point-count workflow searches
% the mode-root grids returned by `pointsFromModeRoots` and selects the
% largest candidate modal prefix whose grid has exactly the requested
% number of points. `nPoints` and `z` are mutually exclusive.
%
% Unless `weights` are supplied with `z`, the quadrature weights are fitted
% once for the complete candidate band by `quadratureWeightsForPoints`.
% Every leading modal prefix is then assessed using those same physical
% points and weights. The normalized-Gram policy is always active. Optional
% rejected-mode leakage and scalar quadratic-aliasing policies are enabled
% by supplying their positive tolerances. The returned production transform
% contains the largest consecutive prefix accepted by every enabled policy.
%
% If `nModes` is supplied with `z`, that band is strict: every requested
% mode must pass every enabled policy, or construction throws an error.
% With `nModes` omitted, the candidate band is the largest available prefix
% that can be represented by the supplied points and may be reduced by the
% policies. The `nPoints` workflow determines its candidate band from the
% exact mode-root grid and therefore does not accept `nModes`.
%
% For prefix $$N$$, define
%
% $$
% E_N=S_N(\Gamma_N-\Gamma_{0,N})S_N,\qquad
% S_N=\operatorname{diag}\!\left(
% \left|\operatorname{diag}\Gamma_{0,N}\right|^{-1/2}
% \right).
% $$
%
% The always-enabled Gram policy accepts
% $$\lVert E_N\rVert_2\leq\texttt{gramTolerance}$$. The shipped default
% is $$10^{-2}$$. Regression sweeps over constant and representative
% exponential stratifications in physical, WKB, and density coordinates
% found constant cosine rules near roundoff and ordinary exponential rules
% below approximately $$5\times10^{-3}$$; the factor-two margin preserves
% those bands while still rejecting the order-one DCT-I Nyquist failure.
%
% `leakageTolerance` enables rejected-mode leakage
%
% $$
% \ell_N=\max_{N<j\leq N_\mathrm{check}}
% \frac{\lVert\Pi_N^\mathrm{discrete}u_j\rVert_\mu}
% {\lVert u_j\rVert_\mu}.
% $$
%
% Its default check count is twice the candidate count. An explicit
% `nCheckModes` must exceed the candidate band so every prefix has at least
% one rejected comparison mode. `quadraticAliasingTolerance` enables
%
% $$
% q_N=\max_{1\leq i\leq j\leq N}
% \frac{\left\lVert\Pi_N^\mathrm{discrete}(u_i u_j)-
% \Pi_N^\mathrm{continuous}(u_i u_j)\right\rVert_\mu}
% {\lVert u_i u_j\rVert_\mu}.
% $$
%
% The continuous quadratic projection is evaluated on the source solver's
% inner-product grid, independently of the fixed transform rule. Leakage
% and quadratic aliasing require a positive-definite target Gram matrix;
% signed targets support Gram assessment only. `IMDiscreteTransformAssessment`
% defines every table column and policy-result field.
%
% ```matlab
% [transform,assessment] = basisSet.discreteTransform(nPoints=32);
% [transform,assessment] = basisSet.discreteTransform(z=z,nModes=8);
% [transform,assessment] = basisSet.discreteTransform(z=z,weights=w,nModes=8);
% assessment.prefixDiagnostics
% ```
%
% - Topic: Build discrete transforms
% - Declaration: [transform,assessment] = discreteTransform(basisSet,options)
% - Parameter options.nPoints: exact requested point count for a mode-root grid
% - Parameter options.z: increasing explicit physical sample points
% - Parameter options.weights: optional fixed quadrature weights aligned with `z`
% - Parameter options.nModes: optional strict number of leading candidate modes for explicit `z`
% - Parameter options.gramTolerance: nonnegative normalized-Gram operator-error tolerance
% - Parameter options.leakageTolerance: optional positive rejected-mode leakage tolerance
% - Parameter options.quadraticAliasingTolerance: optional positive scalar quadratic-aliasing tolerance
% - Parameter options.nCheckModes: optional number of source modes used by the leakage policy
% - Returns transform: retained production transform
% - Returns assessment: fixed-rule candidate and retained-prefix diagnostics
arguments
    self IMBasisSet
    options.nPoints double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.z (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.weights (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nModes double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.gramTolerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-2
    options.leakageTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.quadraticAliasingTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nCheckModes double {mustBeReal, mustBeFinite} = zeros(0,1)
end

validateOptionalPositiveInteger(options.nPoints,"nPoints");
validateOptionalPositiveInteger(options.nModes,"nModes");
validateOptionalPositiveInteger(options.nCheckModes,"nCheckModes");
validateOptionalPositiveTolerance(options.leakageTolerance,"leakageTolerance");
validateOptionalPositiveTolerance(options.quadraticAliasingTolerance,"quadraticAliasingTolerance");

hasPointCount = ~isempty(options.nPoints);
hasExplicitPoints = ~isempty(options.z);
hasExplicitWeights = ~isempty(options.weights);
hasExplicitModeCount = ~isempty(options.nModes);
if hasPointCount == hasExplicitPoints
    error("IMBasisSet:InvalidDiscretePointSpecification", "Specify exactly one of nPoints or z.");
end
if hasExplicitWeights && ~hasExplicitPoints
    error("IMBasisSet:InvalidDiscretePointSpecification", "weights can be supplied only with explicit z.");
end
if hasPointCount && (hasExplicitWeights || hasExplicitModeCount)
    error("IMBasisSet:InvalidDiscretePointSpecification", "The nPoints workflow determines its mode-root candidate band and cannot be combined with weights or nModes.");
end

nAvailableModes = size(self.nativeModes,2);
if hasPointCount
    requestedPointCount = options.nPoints;
    [z,candidateModeCount] = pointsForExactCount(self,requestedPointCount,nAvailableModes);
else
    z = options.z(:);
    requestedPointCount = length(z);
    if hasExplicitModeCount
        candidateModeCount = options.nModes;
    else
        candidateModeCount = min(nAvailableModes,length(z));
    end
end
if candidateModeCount > nAvailableModes
    error("IMBasisSet:InvalidDiscreteModeCount", "The basis set contains %d modes, but nModes=%d was requested.", nAvailableModes, candidateModeCount);
end

if hasExplicitWeights
    weightFit = [];
    candidateTransform = self.buildDiscreteTransform(z,options.weights(:),candidateModeCount);
else
    [~,weightFit] = self.quadratureWeightsForPoints(z=z,nModes=candidateModeCount);
    candidateTransform = weightFit.transform;
end
if (~isempty(options.leakageTolerance) || ~isempty(options.quadraticAliasingTolerance)) ...
        && ~candidateTransform.targetGramIsPositiveDefinite
    error("IMBasisSet:UnavailableDiscreteTransformPolicy", ...
        "Rejected-mode leakage and quadratic-aliasing policies require a positive-definite target Gram matrix; use the Gram policy alone for signed or indefinite targets.");
end

prefixTransforms = cell(candidateModeCount,1);
gramError = zeros(candidateModeCount,1);
roundTripError = zeros(candidateModeCount,1);
inverseConditionNumber = zeros(candidateModeCount,1);
gramConditionNumber = zeros(candidateModeCount,1);
for nPrefix = 1:candidateModeCount
    prefixTransforms{nPrefix} = prefixTransform(candidateTransform,nPrefix);
    gramError(nPrefix) = prefixTransforms{nPrefix}.relativeGramOperatorError;
    roundTripError(nPrefix) = prefixTransforms{nPrefix}.roundTripError;
    inverseConditionNumber(nPrefix) = prefixTransforms{nPrefix}.inverseMatrixConditionNumber;
    gramConditionNumber(nPrefix) = prefixTransforms{nPrefix}.gramConditionNumber;
end
gramAccepted = cumulativeAcceptance(gramError <= options.gramTolerance);
gramPolicy = policyResult("gram",true,options.gramTolerance,gramError,gramAccepted,candidateModeCount);

leakageError = nan(candidateModeCount,1);
leakageLimitingModeNumber = nan(candidateModeCount,1);
if isempty(options.leakageTolerance)
    leakageAccepted = true(candidateModeCount,1);
    leakagePolicy = policyResult("leakage",false,[],leakageError,leakageAccepted,candidateModeCount);
    leakagePolicy.nCheckModes = [];
    leakagePolicy.limitingModeNumber = leakageLimitingModeNumber;
else
    if isempty(options.nCheckModes)
        nCheckModes = 2*candidateModeCount;
    else
        nCheckModes = options.nCheckModes;
    end
    if nCheckModes <= candidateModeCount
        error("IMBasisSet:InvalidLeakageCheckModeCount", "nCheckModes must be greater than the candidate mode count when leakage assessment is enabled.");
    end
    checkBasis = self;
    if nCheckModes > nAvailableModes
        try
            checkBasis = self.solver.solveEVP(self.evp,nModes=nCheckModes);
        catch cause
            exception = MException("IMBasisSet:AuxiliaryModeUnavailable", "The leakage policy could not obtain %d source modes.", nCheckModes);
            throw(addCause(exception,cause))
        end
        if size(checkBasis.nativeModes,2) < nCheckModes
            error("IMBasisSet:AuxiliaryModeUnavailable", "The leakage policy requested %d source modes, but the solver returned %d.", nCheckModes, size(checkBasis.nativeModes,2));
        end
        referenceEigenvalues = self.eigenvalues(1:candidateModeCount);
        checkEigenvalues = checkBasis.eigenvalues(1:candidateModeCount);
        eigenvalueScale = max(1,max(abs([referenceEigenvalues checkEigenvalues]),[],"all"));
        if ~isequal(checkBasis.modeNumber(1:candidateModeCount),self.modeNumber(1:candidateModeCount)) ...
                || any(abs(checkEigenvalues-referenceEigenvalues) > 1e-8*eigenvalueScale)
            error("IMBasisSet:AuxiliaryModeMismatch", "The leakage-policy solve did not reproduce the candidate modal prefix.");
        end
    end
    checkValues = checkBasis.u(candidateTransform.z);
    checkValues = checkValues(:,1:nCheckModes);
    checkGram = checkBasis.gramMatrix();
    checkNorms = diag(checkGram(1:nCheckModes,1:nCheckModes));
    if any(checkNorms <= 0)
        error("IMBasisSet:UnavailableDiscreteTransformPolicy", ...
            "Rejected-mode leakage requires positive continuous norms for every checked source mode.");
    end
    checkModeNumber = checkBasis.modeNumber(1:nCheckModes);
    [leakageError,leakageLimitingModeNumber] = rejectedModeLeakage(prefixTransforms,checkValues,checkNorms,checkModeNumber);
    leakageAccepted = cumulativeAcceptance(leakageError <= options.leakageTolerance);
    leakagePolicy = policyResult("leakage",true,options.leakageTolerance,leakageError,leakageAccepted,candidateModeCount);
    leakagePolicy.nCheckModes = nCheckModes;
    leakagePolicy.limitingModeNumber = leakageLimitingModeNumber;
end

quadraticError = nan(candidateModeCount,1);
quadraticLimitingModeNumberI = nan(candidateModeCount,1);
quadraticLimitingModeNumberJ = nan(candidateModeCount,1);
if isempty(options.quadraticAliasingTolerance)
    quadraticAccepted = true(candidateModeCount,1);
    quadraticPolicy = policyResult("quadraticAliasing",false,[],quadraticError,quadraticAccepted,candidateModeCount);
    quadraticPolicy.limitingModeNumberI = quadraticLimitingModeNumberI;
    quadraticPolicy.limitingModeNumberJ = quadraticLimitingModeNumberJ;
else
    integrationGrid = self.solver.innerProductGrid(self.zDomain);
    integrationModes = self.u(integrationGrid);
    integrationModes = integrationModes(:,1:candidateModeCount);
    context = self.evp.contextForSolver(self.solver);
    spec = self.evp.innerProduct();
    interiorWeight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,integrationGrid,context);
    if isscalar(interiorWeight)
        interiorWeight = interiorWeight*ones(size(integrationGrid));
    else
        interiorWeight = interiorWeight(:);
    end
    if numel(interiorWeight) ~= length(integrationGrid) || ~isreal(interiorWeight) || any(~isfinite(interiorWeight))
        error("IMBasisSet:InvalidDiscreteMetricWeight", "The EVP interior weight must return one finite real value for each inner-product point.");
    end
    endpointModeValues = zeros(0,candidateModeCount);
    endpointCoefficients = zeros(0,1);
    endpointWeights = [spec.surfaceWeights; spec.bottomWeights];
    for iWeight = 1:numel(endpointWeights)
        endpointWeight = endpointWeights(iWeight);
        if endpointWeight.d ~= 0
            error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "Scalar quadratic-aliasing assessment requires endpoint terms that depend only on values.");
        end
        if endpointWeight.c == 0
            continue
        end
        if string(endpointWeight.location) == "surface"
            zEndpoint = self.zDomain(2);
        else
            zEndpoint = self.zDomain(1);
        end
        valuesAtEndpoint = self.u(zEndpoint);
        endpointModeValues(end+1,:) = valuesAtEndpoint(1,1:candidateModeCount); %#ok<AGROW>
        endpointCoefficients(end+1,1) = endpointWeight.coefficient*endpointWeight.c*endpointWeight.c; %#ok<AGROW>
    end
    integrate = @(integrand) self.solver.integrateInnerProduct(integrationGrid,integrand,self.zDomain);
    [quadraticError,quadraticLimitingModeNumberI,quadraticLimitingModeNumberJ] = scalarQuadraticAliasing( ...
        candidateTransform,prefixTransforms,integrationModes,interiorWeight,integrate,endpointModeValues,endpointCoefficients);
    quadraticAccepted = cumulativeAcceptance(quadraticError <= options.quadraticAliasingTolerance);
    quadraticPolicy = policyResult("quadraticAliasing",true,options.quadraticAliasingTolerance,quadraticError,quadraticAccepted,candidateModeCount);
    quadraticPolicy.limitingModeNumberI = quadraticLimitingModeNumberI;
    quadraticPolicy.limitingModeNumberJ = quadraticLimitingModeNumberJ;
end

combinedAccepted = gramAccepted & leakageAccepted & quadraticAccepted;
retainedModeCount = find(combinedAccepted,1,"last");
if isempty(retainedModeCount)
    if hasExplicitModeCount
        throwStrictModeFailure(candidateModeCount,gramPolicy,leakagePolicy,quadraticPolicy);
    end
    error("IMBasisSet:NoAcceptableDiscreteTransformPrefix", "No candidate mode passes every enabled retained-band policy on this fixed quadrature rule.");
end
if hasExplicitModeCount && retainedModeCount < candidateModeCount
    throwStrictModeFailure(candidateModeCount,gramPolicy,leakagePolicy,quadraticPolicy);
end

[limitingPolicy,retentionReason] = combinedPolicyResult(candidateModeCount,retainedModeCount,gramPolicy,leakagePolicy,quadraticPolicy);
transform = prefixTransforms{retainedModeCount};
modeCount = (1:candidateModeCount).';
lastModeNumber = candidateTransform.modeNumber(:);
prefixDiagnostics = table(modeCount,lastModeNumber,gramError,roundTripError,inverseConditionNumber,gramConditionNumber, ...
    leakageError,leakageLimitingModeNumber,quadraticError,quadraticLimitingModeNumberI,quadraticLimitingModeNumberJ, ...
    gramAccepted,leakageAccepted,quadraticAccepted,combinedAccepted, ...
    VariableNames=["modeCount","lastModeNumber","gramError","roundTripError","inverseMatrixConditionNumber","gramConditionNumber", ...
    "leakageError","leakageLimitingModeNumber","quadraticAliasingError","quadraticLimitingModeNumberI","quadraticLimitingModeNumberJ", ...
    "gramAccepted","leakageAccepted","quadraticAccepted","combinedAccepted"]);
assessment = IMDiscreteTransformAssessment(transform=transform,candidateTransform=candidateTransform,weightFit=weightFit, ...
    requestedPointCount=requestedPointCount,prefixDiagnostics=prefixDiagnostics,gramPolicy=gramPolicy, ...
    leakagePolicy=leakagePolicy,quadraticAliasingPolicy=quadraticPolicy,limitingPolicy=limitingPolicy,retentionReason=retentionReason);
end

function validateOptionalPositiveInteger(value,name)
if isempty(value)
    return
end
if ~isscalar(value) || value <= 0 || value ~= floor(value)
    error("IMBasisSet:InvalidDiscreteTransformOption", "%s must be empty or a positive integer scalar.", name);
end
end

function validateOptionalPositiveTolerance(value,name)
if isempty(value)
    return
end
if ~isscalar(value) || value <= 0
    error("IMBasisSet:InvalidDiscreteTransformOption", "%s must be empty or a positive finite scalar.", name);
end
end

function [z,candidateModeCount] = pointsForExactCount(basisSet,nPoints,nAvailableModes)
pointCounts = nan(nAvailableModes,1);
pointGrids = cell(nAvailableModes,1);
maximumPossibleModeCount = min(nAvailableModes,nPoints);
for nModes = maximumPossibleModeCount:-1:1
    try
        pointGrids{nModes} = basisSet.pointsFromModeRoots(nModes=nModes);
        pointCounts(nModes) = length(pointGrids{nModes});
    catch exception
        if nModes == nAvailableModes && strcmp(exception.identifier,"IMBasisSet:AuxiliaryModeUnavailable")
            continue
        end
        rethrow(exception)
    end
    if pointCounts(nModes) == nPoints
        candidateModeCount = nModes;
        z = pointGrids{nModes};
        return
    end
end
for nModes = (maximumPossibleModeCount+1):nAvailableModes
    try
        pointGrids{nModes} = basisSet.pointsFromModeRoots(nModes=nModes);
        pointCounts(nModes) = length(pointGrids{nModes});
    catch exception
        if nModes == nAvailableModes && strcmp(exception.identifier,"IMBasisSet:AuxiliaryModeUnavailable")
            continue
        end
        rethrow(exception)
    end
end
attainableCounts = unique(pointCounts(isfinite(pointCounts)));
[~,order] = sort(abs(attainableCounts-nPoints));
nearestCounts = attainableCounts(order(1:min(3,length(order))));
nearestText = join(string(nearestCounts(:).'),", ");
error("IMBasisSet:UnattainableDiscretePointCount", ...
    "No available mode-root grid has exactly %d points. Nearest attainable counts are %s.",nPoints,nearestText);
end

function transform = prefixTransform(candidateTransform,nModes)
transform = IMDiscreteTransform(z=candidateTransform.z,weights=candidateTransform.weights, ...
    modeNumber=candidateTransform.modeNumber(1:nModes),normalization=candidateTransform.normalization, ...
    inverseMatrix=candidateTransform.inverseMatrix(:,1:nModes),metricMatrix=candidateTransform.metricMatrix, ...
    targetGramMatrix=candidateTransform.targetGramMatrix(1:nModes,1:nModes));
end

function accepted = cumulativeAcceptance(rawAccepted)
accepted = cumprod(double(rawAccepted(:))) > 0;
end

function policy = policyResult(name,enabled,tolerance,errorValues,accepted,candidateModeCount)
acceptedModeCount = find(accepted,1,"last");
if isempty(acceptedModeCount)
    acceptedModeCount = 0;
end
if ~enabled
    limitingValue = NaN;
    reason = sprintf("The %s policy is disabled.", name);
elseif acceptedModeCount < candidateModeCount
    limitingValue = errorValues(acceptedModeCount+1);
    reason = sprintf("Accepted %d of %d candidate modes; the next prefix has value %.6g above tolerance %.6g.", ...
        acceptedModeCount,candidateModeCount,limitingValue,tolerance);
else
    limitingValue = max(errorValues);
    reason = sprintf("All %d candidate modes pass tolerance %.6g; the largest policy value is %.6g.", ...
        candidateModeCount,tolerance,limitingValue);
end
policy = struct(name=string(name),enabled=logical(enabled),tolerance=tolerance,error=errorValues(:), ...
    accepted=accepted(:),acceptedModeCount=acceptedModeCount,maximumAcceptedModeCount=acceptedModeCount, ...
    limitingValue=limitingValue,reason=string(reason),limited=enabled && acceptedModeCount < candidateModeCount);
end

function [errors,limitingModeNumber] = rejectedModeLeakage(prefixTransforms,checkValues,checkNorms,checkModeNumber)
nCandidate = length(prefixTransforms);
nCheckModes = length(checkModeNumber);
errors = zeros(nCandidate,1);
limitingModeNumber = nan(nCandidate,1);
for nPrefix = 1:nCandidate
    transform = prefixTransforms{nPrefix};
    rejectedIndices = (nPrefix+1):nCheckModes;
    coefficients = transform.forwardMatrix*checkValues(:,rejectedIndices);
    targetNorms = abs(diag(transform.targetGramMatrix));
    numeratorSquared = sum(targetNorms.*abs(coefficients).^2,1);
    denominatorSquared = abs(checkNorms(rejectedIndices)).';
    leakage = sqrt(numeratorSquared./denominatorSquared);
    leakage(~isfinite(leakage)) = Inf;
    [errors(nPrefix),iLimiting] = max(leakage);
    limitingModeNumber(nPrefix) = checkModeNumber(rejectedIndices(iLimiting));
end
end

function [errors,limitingModeNumberI,limitingModeNumberJ] = scalarQuadraticAliasing(candidateTransform,prefixTransforms, ...
        integrationModes,interiorWeight,integrate,endpointModeValues,endpointCoefficients)
nCandidate = length(prefixTransforms);
nPairs = nCandidate*(nCandidate+1)/2;
pairIndices = zeros(nPairs,2);
sampleProducts = zeros(length(candidateTransform.z),nPairs);
continuousPairings = zeros(nCandidate,nPairs);
productNorms = zeros(nPairs,1);
iPair = 0;
for iMode = 1:nCandidate
    for jMode = iMode:nCandidate
        iPair = iPair + 1;
        pairIndices(iPair,:) = [iMode jMode];
        sampleProducts(:,iPair) = candidateTransform.inverseMatrix(:,iMode).*candidateTransform.inverseMatrix(:,jMode);
        product = integrationModes(:,iMode).*integrationModes(:,jMode);
        productNorm = integrate(interiorWeight.*product.*product);
        for iEndpoint = 1:length(endpointCoefficients)
            productAtEndpoint = endpointModeValues(iEndpoint,iMode)*endpointModeValues(iEndpoint,jMode);
            productNorm = productNorm + endpointCoefficients(iEndpoint)*productAtEndpoint*productAtEndpoint;
        end
        if productNorm <= 100*eps(max(1,abs(productNorm)))
            error("IMBasisSet:UnavailableDiscreteTransformPolicy", ...
                "Scalar quadratic-aliasing assessment requires a positive continuous norm for every retained-mode product.");
        end
        productNorms(iPair) = productNorm;
        for kMode = 1:nCandidate
            pairing = integrate(interiorWeight.*integrationModes(:,kMode).*product);
            for iEndpoint = 1:length(endpointCoefficients)
                productAtEndpoint = endpointModeValues(iEndpoint,iMode)*endpointModeValues(iEndpoint,jMode);
                pairing = pairing + endpointCoefficients(iEndpoint)*endpointModeValues(iEndpoint,kMode)*productAtEndpoint;
            end
            continuousPairings(kMode,iPair) = pairing;
        end
    end
end

errors = zeros(nCandidate,1);
limitingModeNumberI = nan(nCandidate,1);
limitingModeNumberJ = nan(nCandidate,1);
for nPrefix = 1:nCandidate
    transform = prefixTransforms{nPrefix};
    includedPairs = find(pairIndices(:,1) <= nPrefix & pairIndices(:,2) <= nPrefix);
    targetNorms = diag(transform.targetGramMatrix);
    discreteCoefficients = transform.forwardMatrix*sampleProducts(:,includedPairs);
    continuousCoefficients = continuousPairings(1:nPrefix,includedPairs)./targetNorms;
    coefficientDifference = discreteCoefficients-continuousCoefficients;
    numeratorSquared = sum(abs(targetNorms).*abs(coefficientDifference).^2,1);
    denominatorSquared = productNorms(includedPairs).';
    pairErrors = sqrt(numeratorSquared./denominatorSquared);
    pairErrors(~isfinite(pairErrors)) = Inf;
    [errors(nPrefix),iLimiting] = max(pairErrors);
    limitingPair = pairIndices(includedPairs(iLimiting),:);
    limitingModeNumberI(nPrefix) = candidateTransform.modeNumber(limitingPair(1));
    limitingModeNumberJ(nPrefix) = candidateTransform.modeNumber(limitingPair(2));
end
end

function throwStrictModeFailure(requestedModeCount,gramPolicy,leakagePolicy,quadraticPolicy)
policies = {gramPolicy,leakagePolicy,quadraticPolicy};
failureDetails = strings(0,1);
for iPolicy = 1:length(policies)
    policy = policies{iPolicy};
    if policy.enabled && policy.acceptedModeCount < requestedModeCount
        failureDetails(end+1,1) = sprintf("%s (accepted %d, first failing value %.6g, tolerance %.6g)", ...
            policy.name,policy.acceptedModeCount,policy.limitingValue,policy.tolerance); %#ok<AGROW>
    end
end
error("IMBasisSet:DiscreteTransformPolicyFailed", ...
    "The explicitly requested %d-mode band failed retained-band policy: %s.", requestedModeCount, join(failureDetails,"; "));
end

function [limitingPolicy,reason] = combinedPolicyResult(candidateModeCount,retainedModeCount,gramPolicy,leakagePolicy,quadraticPolicy)
if retainedModeCount == candidateModeCount
    limitingPolicy = "none";
    reason = sprintf("All %d candidate modes pass every enabled retained-band policy.", candidateModeCount);
    return
end
policies = {gramPolicy,leakagePolicy,quadraticPolicy};
limitingNames = strings(0,1);
for iPolicy = 1:length(policies)
    policy = policies{iPolicy};
    if policy.enabled && policy.acceptedModeCount == retainedModeCount
        limitingNames(end+1,1) = policy.name; %#ok<AGROW>
    end
end
limitingPolicy = join(limitingNames," and ");
reason = sprintf("Retained %d of %d candidate modes; the %s policy rejected the next prefix.", ...
    retainedModeCount,candidateModeCount,limitingPolicy);
end
