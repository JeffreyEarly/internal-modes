function [weights, weightFit] = quadratureWeightsForPoints(self, options)
% Find quadrature weights for fixed physical sample points.
%
% The points `z` are fixed. This method solves for one weight $$w_k$$ per
% point so the sampled inner products of the first `nModes` retained modes
% approximate their continuous Gram matrix. Let $$\Gamma(w)$$ be the
% sampled Gram matrix, let $$\Gamma_0$$ be its continuous target, and define
%
% $$
% S=\operatorname{diag}\!\left(
% \left|\operatorname{diag}\Gamma_0\right|^{-1/2}
% \right),
% \qquad
% E(w)=S\left(\Gamma(w)-\Gamma_0\right)S.
% $$
%
% The default `"normalizedGramFrobenius"` objective minimizes
% $$\|E(w)\|_{\mathrm F}$$. It combines the mismatch from every retained
% mode pair: diagonal entries of $$E$$ measure modal norm errors, while
% off-diagonal entries measure lost orthogonality. For target modal norms
% $$C_i=(\Gamma_0)_{ii}$$ and fixed endpoint contribution
% $$\Gamma_{\mathrm{endpoint}}$$, the corresponding least-squares system is
%
% $$
% (A_{\mathrm{LS}})_{(i,j),k}
% =\frac{r(z_k)\Phi_{ki}\Phi_{kj}}{\sqrt{|C_iC_j|}},
% \qquad
% (b_{\mathrm{LS}})_{(i,j)}
% =\frac{(\Gamma_0-\Gamma_{\mathrm{endpoint}})_{ij}}
% {\sqrt{|C_iC_j|}}.
% $$
%
% The fitted weights solve
%
% $$
% \min_w\left\|A_{\mathrm{LS}}w-b_{\mathrm{LS}}\right\|_2.
% $$
%
% By default they also satisfy
%
% $$
% w_k\geq0,
% \qquad
% \sum_k w_k=z_\mathrm{surface}-z_\mathrm{bottom}.
% $$
%
% Geometric control-volume weights provide the initial guess and reference
% comparison. Unlike those geometric widths, fitted weights are algebraic
% quadrature coefficients and need not remain positive or sum to the full
% depth when the corresponding constraints are disabled. The optimizer
% works with dimensionless weights $$x_k=w_k/D$$, where $$D$$ is the full
% depth. Custom objectives continue to define $$Aw-b$$ in physical units.
% This method requires `lsqlin` from Optimization Toolbox.
%
% For the built-in objective, `weightFit.residualNorm` is the aggregate
% error $$\|E(w)\|_{\mathrm F}$$. The resulting transform separately reports
% `relativeGramOperatorError` as $$\|E(w)\|_2$$, the largest Gram distortion
% over any normalized combination of retained modes. `roundTripError`
% measures algebraic coefficient recovery and can be tiny even when either
% Gram error is appreciable. The quadrature-weight regression sweep supports
% retaining the unregularized Frobenius objective; geometric weights remain
% the initial guess and comparison baseline.
%
% A custom objective is a function handle accepting a context struct and
% returning a scalar struct with fields `A`, `b`, and optional `name`. The
% context contains `z`, `modeNumber`, `normalization`, `inverseMatrix`,
% `interiorWeight`, `targetGramMatrix`, `endpointGramMatrix`,
% `geometricWeights`, `normalizedGramA`, and `normalizedGramB`.
%
% ```matlab
% [weights,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
% transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
% [weightFit.residualNorm weightFit.geometricResidualNorm]
% [transform.relativeGramOperatorError weightFit.geometricTransform.relativeGramOperatorError]
% ```
%
% - Topic: Build discrete transforms
% - Declaration: [weights,weightFit] = quadratureWeightsForPoints(basisSet,options)
% - Parameter options.z: increasing fixed physical sample points
% - Parameter options.nModes: number of leading retained modes
% - Parameter options.objective: `"normalizedGramFrobenius"` or a custom least-squares callback
% - Parameter options.nonnegative: whether fitted weights must be nonnegative
% - Parameter options.constrainDepth: whether fitted weights must sum to the full depth
% - Returns weights: fitted quadrature weights aligned with `z`
% - Returns weightFit: optional fitted-versus-geometric diagnostics
arguments
    self IMBasisSet
    options.z (:,1) double {mustBeReal, mustBeFinite}
    options.nModes (1,1) double {mustBeInteger, mustBePositive} = size(self.nativeModes,2)
    options.objective = "normalizedGramFrobenius"
    options.nonnegative (1,1) logical = true
    options.constrainDepth (1,1) logical = true
end

z = options.z(:);
nModes = options.nModes;
depth = diff(self.zDomain);
edges = [self.zDomain(1); 0.5*(z(1:end-1) + z(2:end)); self.zDomain(2)];
geometricWeights = diff(edges);
geometricTransform = self.discreteTransform(z=z,weights=geometricWeights,nModes=nModes);

evpContext = self.evp.contextForSolver(self.solver);
innerProduct = self.evp.innerProduct();
interiorWeight = IMEigenvalueProblem.evaluateCoefficient(innerProduct.interiorWeight,z,evpContext);
if isscalar(interiorWeight)
    interiorWeight = interiorWeight*ones(size(z));
else
    interiorWeight = interiorWeight(:);
end

inverseMatrix = geometricTransform.inverseMatrix;
targetGramMatrix = geometricTransform.targetGramMatrix;
interiorMetricMatrix = diag(interiorWeight.*geometricWeights);
endpointMetricMatrix = geometricTransform.metricMatrix - interiorMetricMatrix;
endpointGramMatrix = inverseMatrix.'*endpointMetricMatrix*inverseMatrix;
endpointGramMatrix = 0.5*(endpointGramMatrix + endpointGramMatrix.');
targetNorms = diag(targetGramMatrix);

nSamples = length(z);
normalizedGramA = zeros(nModes*nModes,nSamples);
normalizedGramB = zeros(nModes*nModes,1);
iRow = 0;
for iMode = 1:nModes
    for jMode = 1:nModes
        iRow = iRow + 1;
        scale = sqrt(abs(targetNorms(iMode)*targetNorms(jMode)));
        normalizedGramA(iRow,:) = (interiorWeight.*inverseMatrix(:,iMode).*inverseMatrix(:,jMode)).'/scale;
        normalizedGramB(iRow) = (targetGramMatrix(iMode,jMode) - endpointGramMatrix(iMode,jMode))/scale;
    end
end

objectiveContext = struct();
objectiveContext.z = z;
objectiveContext.modeNumber = geometricTransform.modeNumber;
objectiveContext.normalization = geometricTransform.normalization;
objectiveContext.inverseMatrix = inverseMatrix;
objectiveContext.interiorWeight = interiorWeight;
objectiveContext.targetGramMatrix = targetGramMatrix;
objectiveContext.endpointGramMatrix = endpointGramMatrix;
objectiveContext.geometricWeights = geometricWeights;
objectiveContext.normalizedGramA = normalizedGramA;
objectiveContext.normalizedGramB = normalizedGramB;

objective = options.objective;
if ischar(objective) || (isstring(objective) && isscalar(objective))
    objectiveName = string(objective);
    if objectiveName ~= "normalizedGramFrobenius"
        error("IMBasisSet:UnknownQuadratureObjective", "Unknown quadrature objective ""%s"".", objectiveName);
    end
    objectiveMatrix = normalizedGramA;
    objectiveTarget = normalizedGramB;
elseif isa(objective,"function_handle")
    try
        specification = objective(objectiveContext);
    catch cause
        exception = MException("IMBasisSet:QuadratureObjectiveFailed", "The custom quadrature objective failed while building its least-squares system.");
        throw(addCause(exception,cause))
    end
    if ~isstruct(specification) || ~isscalar(specification) || ~isfield(specification,"A") || ~isfield(specification,"b")
        error("IMBasisSet:InvalidQuadratureObjective", "A custom quadrature objective must return a scalar struct with fields A and b.");
    end
    objectiveMatrix = specification.A;
    objectiveTarget = specification.b;
    objectiveName = "custom";
    if isfield(specification,"name")
        if ~(ischar(specification.name) || (isstring(specification.name) && isscalar(specification.name)))
            error("IMBasisSet:InvalidQuadratureObjective", "The optional custom objective name must be a text scalar.");
        end
        objectiveName = string(specification.name);
    end
else
    error("IMBasisSet:InvalidQuadratureObjective", "objective must be ""normalizedGramFrobenius"" or a function handle returning a least-squares system.");
end

if ~isnumeric(objectiveMatrix) || ~ismatrix(objectiveMatrix) || ~isreal(objectiveMatrix) || isempty(objectiveMatrix) ...
        || size(objectiveMatrix,2) ~= nSamples || any(~isfinite(objectiveMatrix),"all")
    error("IMBasisSet:InvalidQuadratureObjective", "The custom objective matrix A must be a finite real matrix with one column per sample point.");
end
if ~isnumeric(objectiveTarget) || ~isvector(objectiveTarget) || ~isreal(objectiveTarget) ...
        || length(objectiveTarget) ~= size(objectiveMatrix,1) || any(~isfinite(objectiveTarget),"all")
    error("IMBasisSet:InvalidQuadratureObjective", "The custom objective target b must be a finite real vector with one entry per row of A.");
end
objectiveMatrix = double(objectiveMatrix);
objectiveTarget = double(objectiveTarget(:));

if isempty(which("lsqlin"))
    error("IMBasisSet:MissingQuadratureOptimizer", "Fitting quadrature weights requires lsqlin from Optimization Toolbox. Supply weights explicitly when Optimization Toolbox is unavailable.");
end

if options.nonnegative
    lowerBounds = zeros(nSamples,1);
else
    lowerBounds = [];
end
if options.constrainDepth
    equalityMatrix = ones(1,nSamples);
    equalityTarget = 1;
else
    equalityMatrix = [];
    equalityTarget = [];
end

scaledObjectiveMatrix = depth*objectiveMatrix;
initialScaledWeights = geometricWeights/depth;
solverOptions = optimoptions("lsqlin",Algorithm="active-set",Display="off",ConstraintTolerance=1e-12,OptimalityTolerance=1e-12,StepTolerance=1e-14,MaxIterations=1000);
try
    [scaledWeights,~,~,exitFlag,solverOutput] = lsqlin(scaledObjectiveMatrix,objectiveTarget,[],[],equalityMatrix,equalityTarget,lowerBounds,[],initialScaledWeights,solverOptions);
catch cause
    exception = MException("IMBasisSet:QuadratureWeightFitFailed", "lsqlin failed while fitting quadrature weights on the supplied points.");
    throw(addCause(exception,cause))
end
if isempty(scaledWeights) || exitFlag <= 0 || any(~isfinite(scaledWeights))
    error("IMBasisSet:QuadratureWeightFitFailed", "lsqlin did not converge to finite quadrature weights (exit flag %d).", exitFlag);
end
weights = depth*scaledWeights;

constraintTolerance = 1e-10*max(1,depth);
if options.nonnegative && any(weights < -constraintTolerance)
    error("IMBasisSet:QuadratureWeightFitFailed", "The fitted weights violate the requested nonnegative constraint.");
end
if options.constrainDepth && abs(sum(weights) - depth) > constraintTolerance
    error("IMBasisSet:QuadratureWeightFitFailed", "The fitted weights violate the requested full-depth constraint.");
end
weights(weights < 0 & weights >= -constraintTolerance) = 0;

if nargout > 1
    transform = self.discreteTransform(z=z,weights=weights,nModes=nModes);
    weightFit = IMQuadratureWeightFit(transform=transform,geometricTransform=geometricTransform,objectiveName=objectiveName,objectiveMatrix=objectiveMatrix,objectiveTarget=objectiveTarget,nonnegativeConstraint=options.nonnegative,depthConstraint=options.constrainDepth,depthTarget=depth,exitFlag=exitFlag,solverOutput=solverOutput);
end
end
