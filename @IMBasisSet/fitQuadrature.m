function fit = fitQuadrature(self, options)
% Fit quadrature increments on fixed physical sample points.
%
% The default objective fits the normalized scalar Gram matrix. For target
% norms $$C_i=(\Gamma_0)_{ii}$$ and a fixed endpoint contribution
% $$\Gamma_{\mathrm{end}}$$, the ordered mode-pair residual system is
%
% $$
% A_{(i,j),k}=\frac{r(z_k)\Phi_{ki}\Phi_{kj}}{\sqrt{|C_iC_j|}},
% \qquad
% b_{(i,j)}=\frac{(\Gamma_0-\Gamma_{\mathrm{end}})_{ij}}
% {\sqrt{|C_iC_j|}}.
% $$
%
% The fitted increments minimize $$\|A\Delta z-b\|_2$$. By default they
% also satisfy $$\Delta z_k\geq0$$ and
% $$\sum_k\Delta z_k=z_\mathrm{surface}-z_\mathrm{bottom}$$. The result
% compares the fitted transform with geometric control-volume increments.
% This method requires `lsqlin` from Optimization Toolbox.
%
% A custom objective is a function handle accepting a context struct and
% returning a scalar struct with fields `A`, `b`, and optional `name`. The
% context contains `z`, `modeNumber`, `normalization`, `basisMatrix`,
% `interiorWeight`, `targetGramMatrix`, `endpointGramMatrix`,
% `geometricIncrements`, `normalizedGramA`, and `normalizedGramB`.
%
% ```matlab
% fit = basisSet.fitQuadrature(z=z,nModes=8);
% transform = fit.fittedTransform;
% ```
%
% - Topic: Build discrete transforms
% - Declaration: fit = fitQuadrature(basisSet,options)
% - Parameter options.z: increasing fixed physical sample points
% - Parameter options.nModes: number of leading retained modes
% - Parameter options.objective: `"normalizedGram"` or a custom least-squares callback
% - Parameter options.nonnegative: whether fitted increments must be nonnegative
% - Parameter options.constrainDepth: whether fitted increments must sum to the full depth
% - Returns fit: fitted and geometric quadrature comparison
arguments
    self IMBasisSet
    options.z (:,1) double {mustBeReal, mustBeFinite}
    options.nModes (1,1) double {mustBeInteger, mustBePositive} = size(self.nativeModes,2)
    options.objective = "normalizedGram"
    options.nonnegative (1,1) logical = true
    options.constrainDepth (1,1) logical = true
end

z = options.z(:);
nModes = options.nModes;
depth = diff(self.zDomain);
edges = [self.zDomain(1); 0.5*(z(1:end-1) + z(2:end)); self.zDomain(2)];
geometricIncrements = diff(edges);
geometricTransform = self.discreteTransform(z=z, increments=geometricIncrements, nModes=nModes);

evpContext = self.evp.contextForSolver(self.solver);
innerProduct = self.evp.innerProduct();
interiorWeight = IMEigenvalueProblem.evaluateCoefficient(innerProduct.interiorWeight, z, evpContext);
if isscalar(interiorWeight)
    interiorWeight = interiorWeight*ones(size(z));
else
    interiorWeight = interiorWeight(:);
end

basisMatrix = geometricTransform.basisMatrix;
targetGramMatrix = geometricTransform.targetGramMatrix;
interiorMetricMatrix = diag(interiorWeight.*geometricIncrements);
endpointMetricMatrix = geometricTransform.metricMatrix - interiorMetricMatrix;
endpointGramMatrix = basisMatrix.'*endpointMetricMatrix*basisMatrix;
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
        normalizedGramA(iRow,:) = (interiorWeight.*basisMatrix(:,iMode).*basisMatrix(:,jMode)).'/scale;
        normalizedGramB(iRow) = (targetGramMatrix(iMode,jMode) - endpointGramMatrix(iMode,jMode))/scale;
    end
end

objectiveContext = struct();
objectiveContext.z = z;
objectiveContext.modeNumber = geometricTransform.modeNumber;
objectiveContext.normalization = geometricTransform.normalization;
objectiveContext.basisMatrix = basisMatrix;
objectiveContext.interiorWeight = interiorWeight;
objectiveContext.targetGramMatrix = targetGramMatrix;
objectiveContext.endpointGramMatrix = endpointGramMatrix;
objectiveContext.geometricIncrements = geometricIncrements;
objectiveContext.normalizedGramA = normalizedGramA;
objectiveContext.normalizedGramB = normalizedGramB;

objective = options.objective;
if ischar(objective) || (isstring(objective) && isscalar(objective))
    objectiveName = string(objective);
    if objectiveName ~= "normalizedGram"
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
    error("IMBasisSet:InvalidQuadratureObjective", "objective must be ""normalizedGram"" or a function handle returning a least-squares system.");
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
    error("IMBasisSet:MissingQuadratureOptimizer", "Fitting quadrature increments requires lsqlin from Optimization Toolbox. Supply increments explicitly when Optimization Toolbox is unavailable.");
end

if options.nonnegative
    lowerBounds = zeros(nSamples,1);
else
    lowerBounds = [];
end
if options.constrainDepth
    equalityMatrix = ones(1,nSamples);
    equalityTarget = depth;
else
    equalityMatrix = [];
    equalityTarget = [];
end

solverOptions = optimoptions("lsqlin",Algorithm="active-set",Display="off",ConstraintTolerance=1e-12, ...
    OptimalityTolerance=1e-12,StepTolerance=1e-14,MaxIterations=1000);
try
    [fittedIncrements,~,~,exitFlag,solverOutput] = lsqlin(objectiveMatrix, objectiveTarget, [], [], equalityMatrix, equalityTarget, lowerBounds, [], geometricIncrements, solverOptions);
catch cause
    exception = MException("IMBasisSet:QuadratureFitFailed", "lsqlin failed while fitting quadrature increments on the supplied points.");
    throw(addCause(exception,cause))
end
if isempty(fittedIncrements) || exitFlag <= 0 || any(~isfinite(fittedIncrements))
    error("IMBasisSet:QuadratureFitFailed", "lsqlin did not converge to a finite quadrature fit (exit flag %d).", exitFlag);
end

constraintTolerance = 1e-10*max(1,depth);
if options.nonnegative && any(fittedIncrements < -constraintTolerance)
    error("IMBasisSet:QuadratureFitFailed", "The fitted increments violate the requested nonnegative constraint.");
end
if options.constrainDepth && abs(sum(fittedIncrements) - depth) > constraintTolerance
    error("IMBasisSet:QuadratureFitFailed", "The fitted increments violate the requested full-depth constraint.");
end
fittedIncrements(fittedIncrements < 0 & fittedIncrements >= -constraintTolerance) = 0;

fittedTransform = self.discreteTransform(z=z, increments=fittedIncrements, nModes=nModes);
fit = IMQuadratureFit(fittedTransform=fittedTransform, geometricTransform=geometricTransform, objectiveName=objectiveName, ...
    objectiveMatrix=objectiveMatrix, objectiveTarget=objectiveTarget, nonnegativeConstraint=options.nonnegative, ...
    depthConstraint=options.constrainDepth, depthTarget=depth, exitFlag=exitFlag, solverOutput=solverOutput);
end
