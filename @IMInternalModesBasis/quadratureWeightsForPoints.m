function [weights, weightFit] = quadratureWeightsForPoints(self, options)
% Fit one quadrature rule to aligned internal-mode F/G channels.
%
% Every requested directly representable channel contributes its compressed
% normalized-Gram Frobenius system. The systems are stacked without an
% additional channel weight and solved for one nonnegative, full-depth
% constrained weight vector by default. Identically zero variable columns
% are omitted from that variable's objective while their family positions
% remain present in the returned transform.
%
% A custom objective receives the stacked system in
% `context.normalizedGramA` and `context.normalizedGramB`, row provenance in
% `context.normalizedGramVariables`, and scalar-style per-variable contexts
% in `context.variableContexts.F` and `.G`.
%
% - Topic: Build discrete transforms
% - Declaration: [weights,weightFit] = quadratureWeightsForPoints(basisSet,options)
% - Parameter options.z: increasing fixed physical sample points
% - Parameter options.nModes: number of leading aligned family modes
% - Parameter options.variables: requested direct channels, F and/or G
% - Parameter options.objective: built-in stacked objective or custom callback
% - Parameter options.nonnegative: impose nonnegative shared weights
% - Parameter options.constrainDepth: impose full-depth shared weights
% - Returns weights: fitted shared quadrature weights
% - Returns weightFit: specialized shared-fit diagnostics
arguments
    self IMInternalModesBasis
    options.z (:,1) double {mustBeReal, mustBeFinite}
    options.nModes (1,1) double {mustBeInteger, mustBePositive} = size(self.nativeModes,2)
    options.variables (1,:) string = strings(1,0)
    options.objective = "stackedNormalizedGramFrobenius"
    options.nonnegative (1,1) logical = true
    options.constrainDepth (1,1) logical = true
end

z = options.z(:);
nModes = options.nModes;
depth = diff(self.zDomain);
edges = [self.zDomain(1);0.5*(z(1:end-1)+z(2:end));self.zDomain(2)];
geometricWeights = diff(edges);
if isempty(options.variables)
    variables = self.directlyRepresentableDiscreteVariables(z);
else
    variables = self.canonicalDiscreteVariables(options.variables);
end
if isempty(variables)
    error("IMInternalModesBasis:NoAvailableDiscreteTransformVariable", "Neither F nor G has a directly representable sampled metric on the supplied points.");
end
geometricTransform = self.buildInternalModesDiscreteTransform(z,geometricWeights,nModes,variables);

normalizedGramA = zeros(0,length(z));
normalizedGramB = zeros(0,1);
normalizedGramVariables = strings(0,1);
normalizedGramModePairs = zeros(0,2);
variableContexts = struct();
for variable = variables
    variableContext = geometricTransform.quadratureFitContext(variable=variable);
    variableA = variableContext.normalizedGramA;
    variableB = variableContext.normalizedGramB;
    variablePairs = variableContext.normalizedGramModePairs;
    rows = size(variableA,1);
    normalizedGramA = [normalizedGramA;variableA]; %#ok<AGROW>
    normalizedGramB = [normalizedGramB;variableB]; %#ok<AGROW>
    normalizedGramVariables = [normalizedGramVariables;repmat(variable,rows,1)]; %#ok<AGROW>
    normalizedGramModePairs = [normalizedGramModePairs;variablePairs]; %#ok<AGROW>
    variableContexts.(char(variable)) = variableContext;
end

objectiveContext = struct(z=z,modeNumber=geometricTransform.modeNumber,normalization=geometricTransform.normalization, ...
    geometricWeights=geometricWeights,availableVariables=variables,normalizedGramA=normalizedGramA,normalizedGramB=normalizedGramB, ...
    normalizedGramVariables=normalizedGramVariables,normalizedGramModePairs=normalizedGramModePairs,variableContexts=variableContexts);

objective = options.objective;
objectiveRowVariables = normalizedGramVariables;
objectiveModePairs = normalizedGramModePairs;
if ischar(objective) || (isstring(objective) && isscalar(objective))
    objectiveName = string(objective);
    if ~ismember(objectiveName,["stackedNormalizedGramFrobenius","normalizedGramFrobenius"])
        error("IMBasisSet:UnknownQuadratureObjective", "Unknown internal-mode quadrature objective ""%s"".", objectiveName);
    end
    objectiveName = "stackedNormalizedGramFrobenius";
    objectiveMatrix = normalizedGramA;
    objectiveTarget = normalizedGramB;
elseif isa(objective,"function_handle")
    try
        specification = objective(objectiveContext);
    catch cause
        exception = MException("IMBasisSet:QuadratureObjectiveFailed", "The custom internal-mode quadrature objective failed while building its least-squares system.");
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
    if isfield(specification,"rowVariables")
        objectiveRowVariables = string(specification.rowVariables(:));
    elseif size(objectiveMatrix,1) ~= length(normalizedGramVariables)
        nDefaultRows = length(normalizedGramVariables);
        nCustomRows = size(objectiveMatrix,1) - nDefaultRows;
        if nCustomRows >= 0
            objectiveRowVariables = [normalizedGramVariables;repmat("custom",nCustomRows,1)];
        else
            objectiveRowVariables = repmat("custom",size(objectiveMatrix,1),1);
        end
    end
    if isfield(specification,"modePairs")
        objectiveModePairs = specification.modePairs;
    elseif size(objectiveMatrix,1) ~= size(normalizedGramModePairs,1)
        nDefaultRows = size(normalizedGramModePairs,1);
        nCustomRows = size(objectiveMatrix,1) - nDefaultRows;
        if nCustomRows >= 0
            objectiveModePairs = [normalizedGramModePairs;zeros(nCustomRows,2)];
        else
            objectiveModePairs = zeros(size(objectiveMatrix,1),2);
        end
    end
else
    error("IMBasisSet:InvalidQuadratureObjective", "objective must name the stacked normalized-Gram objective or be a function handle.");
end

nSamples = length(z);
if ~isnumeric(objectiveMatrix) || ~ismatrix(objectiveMatrix) || ~isreal(objectiveMatrix) ...
        || size(objectiveMatrix,2) ~= nSamples || any(~isfinite(objectiveMatrix),"all")
    error("IMBasisSet:InvalidQuadratureObjective", "The objective matrix A must be finite, real, and have one column per sample point.");
end
if ~isnumeric(objectiveTarget) || ~isvector(objectiveTarget) || ~isreal(objectiveTarget) ...
        || length(objectiveTarget) ~= size(objectiveMatrix,1) || any(~isfinite(objectiveTarget),"all")
    error("IMBasisSet:InvalidQuadratureObjective", "The objective target b must be finite, real, and have one entry per objective row.");
end
if length(objectiveRowVariables) ~= size(objectiveMatrix,1) || size(objectiveModePairs,1) ~= size(objectiveMatrix,1)
    error("IMBasisSet:InvalidQuadratureObjective", "Custom objective row provenance must align with its rows.");
end
objectiveMatrix = double(objectiveMatrix);
objectiveTarget = double(objectiveTarget(:));

if ~isempty(objectiveMatrix) && isempty(which("lsqlin"))
    error("IMBasisSet:MissingQuadratureOptimizer", "Fitting quadrature weights requires lsqlin from Optimization Toolbox. Supply weights explicitly when it is unavailable.");
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
if isempty(objectiveMatrix)
    weights = geometricWeights;
    exitFlag = 1;
    solverOutput = struct("message","No active modal Gram rows; returned the geometric rule satisfying the requested constraints.");
else
    scaledObjectiveMatrix = depth*objectiveMatrix;
    initialScaledWeights = geometricWeights/depth;
    solverOptions = sharedQuadratureSolverOptions();
    try
        [scaledWeights,~,~,exitFlag,solverOutput] = lsqlin(scaledObjectiveMatrix,objectiveTarget,[],[],equalityMatrix,equalityTarget,lowerBounds,[],initialScaledWeights,solverOptions);
    catch cause
        exception = MException("IMBasisSet:QuadratureWeightFitFailed", "lsqlin failed while fitting shared internal-mode quadrature weights.");
        throw(addCause(exception,cause))
    end
    if isempty(scaledWeights) || exitFlag <= 0 || any(~isfinite(scaledWeights))
        error("IMBasisSet:QuadratureWeightFitFailed", "lsqlin did not converge to finite quadrature weights (exit flag %d).",exitFlag);
    end
    weights = depth*scaledWeights;
end
constraintTolerance = 1e-10*max(1,depth);
if options.nonnegative && any(weights < -constraintTolerance)
    error("IMBasisSet:QuadratureWeightFitFailed", "The fitted weights violate the requested nonnegative constraint.");
end
if options.constrainDepth && abs(sum(weights)-depth) > constraintTolerance
    error("IMBasisSet:QuadratureWeightFitFailed", "The fitted weights violate the requested full-depth constraint.");
end
weights(weights < 0 & weights >= -constraintTolerance) = 0;

if nargout > 1
    transform = geometricTransform.transformWithWeights(weights);
    weightFit = IMInternalModesQuadratureWeightFit(transform=transform,geometricTransform=geometricTransform,objectiveName=objectiveName, ...
        objectiveMatrix=objectiveMatrix,objectiveTarget=objectiveTarget,objectiveRowVariables=objectiveRowVariables,objectiveModePairs=objectiveModePairs, ...
        nonnegativeConstraint=options.nonnegative,depthConstraint=options.constrainDepth,depthTarget=depth,exitFlag=exitFlag,solverOutput=solverOutput);
end
end

function options = sharedQuadratureSolverOptions()
persistent cachedOptions
if isempty(cachedOptions)
    cachedOptions = optimoptions("lsqlin",Algorithm="active-set",Display="off",ConstraintTolerance=1e-12,OptimalityTolerance=1e-12,StepTolerance=1e-14,MaxIterations=5000);
end
options = cachedOptions;
end
