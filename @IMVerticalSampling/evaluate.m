function [synthesis,record] = evaluate(self,basis,options)
% Evaluate declared vertical observations of an existing continuous basis.
%
% Cell observations use physical dz integration, split at known source and
% caller breakpoints. Successive positive Gauss rules provide an estimated
% error, not a certified bound or an independent qualification of the modes.
% Relative tolerance is measured against the integral of the profile's
% absolute value, rather than a possibly vanishing signed cell integral.
% Cell-average estimates and scales are divided by physical cell width.
%
% Supported derivatives use the basis evaluator's classical convention.
% Interval joins do not add distributional derivative terms. Fixed source
% hints retain the original inventory after selection. Adaptive coordinate
% covariance holds to integration tolerance and floating-point accuracy,
% rather than bitwise equality. Arbitrary intrinsic evaluators may require
% caller-supplied breakpoints or a larger explicitly justified initial order.
%
% No EVP, grid design, weight fitting, physical dual, activity, fitting norm
% or observation/noise norm is prepared here. Pass the returned matrix to
% IMDiscreteTransform.fromObservationMatrices with explicit caller-owned
% coefficient and observation meanings.
%
% - Topic: Evaluate vertical sampling
% - Parameter basis: existing selected/normalized common continuous basis
% - Parameter options.relativeTolerance: estimated-error tolerance relative to each profile's physical L1 integral, default 1e-11
% - Parameter options.absoluteTolerance: nonnegative scalar or one tolerance per row, in synthesis-entry units (observation per coefficient); unit labels do not convert it
% - Parameter options.initialOrder: optional initial Gauss order; otherwise uses immutable original inventory and source resolution
% - Parameter options.maxOrder: maximum Gauss order, default 8192; budget exhaustion throws without returning a matrix
% - Parameter options.breakpoints: additional known physical joins or coefficient discontinuities
% - Returns synthesis: observation-by-column matrix in exact declared order
% - Returns record: observation/basis identities, conventions and estimated integration evidence
arguments (Input)
    self (1,1) IMVerticalSampling
    basis (1,1) IMModeBasis
    options.relativeTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-11
    options.absoluteTolerance (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 0
    options.initialOrder (1,:) double {mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = []
    options.maxOrder (1,1) double {mustBeReal,mustBeFinite,mustBeInteger,mustBePositive} = 8192
    options.breakpoints (1,:) double {mustBeReal,mustBeFinite} = zeros(1,0)
end
arguments (Output)
    synthesis (:,:) double
    record (1,1) struct
end
count = self.observationCount;
columns = basis.nColumns;
observations = self.observations;
if options.relativeTolerance >= 1 || options.maxOrder < 2 || (~isempty(options.initialOrder) && (~isscalar(options.initialOrder) || options.initialOrder < 2))
    error("IMVerticalSampling:InvalidIntegrationOptions","relativeTolerance must be below one and Gauss orders must be scalar integers of at least two.");
end
absoluteTolerance = options.absoluteTolerance;
if isscalar(absoluteTolerance)
    absoluteTolerance = repmat(absoluteTolerance,count,1);
elseif numel(absoluteTolerance) == count
    absoluteTolerance = absoluteTolerance(:);
else
    error("IMVerticalSampling:InvalidIntegrationOptions","Supply a scalar absoluteTolerance or one value per observation row.");
end
if any(observations.zLower < basis.zDomain(1) | observations.zUpper > basis.zDomain(2)) || any(options.breakpoints < basis.zDomain(1) | options.breakpoints > basis.zDomain(2))
    error("IMVerticalSampling:OutsideDomain","Every observation support and caller breakpoint must lie within the physical basis domain.");
end
variables = observations.variable;
orders = observations.derivativeOrder;
for alias = basis.capabilities.aliases
    selected = variables == alias.variable;
    variables(selected) = alias.baseVariable;
    orders(selected) = orders(selected)+alias.derivativeOffset;
end
for row = 1:count
    variableIndex = find(basis.capabilities.variables == variables(row),1);
    if isempty(variableIndex)
        error("IMVerticalSampling:UnavailableVariable","Observation %s requests %s, which the %s family does not provide.",observations.observationLabel(row),variables(row),basis.family);
    end
    if ~ismember(orders(row),basis.capabilities.derivativeOrders{variableIndex})
        error("IMVerticalSampling:UnsupportedDerivative","Observation %s requests unsupported derivative order %d of %s.",observations.observationLabel(row),orders(row),variables(row));
    end
end
info = basis.samplingIntegrationInfo();
breakpoints = unique([info.breakpoints;options.breakpoints(:)]);
breakpoints = breakpoints(breakpoints > basis.zDomain(1) & breakpoints < basis.zDomain(2));
if isempty(options.initialOrder)
    initialOrder = max([128,4*info.originalColumnCount+32,info.sourceResolution]);
    initialOrderPolicy = "max(128,4*originalColumnCount+32,sourceResolution)";
else
    initialOrder = options.initialOrder;
    initialOrderPolicy = "callerProvided";
end
synthesis = zeros(count,columns);
estimatedError = NaN(count,columns);
l1Scale = NaN(count,columns);
finalOrder = zeros(count,1);
subintervalCount = zeros(count,1);
evaluationPointCount = zeros(count,1);
converged = false(count,1);
status = repmat("notApplicable",count,1);
cachedOrders = zeros(1,0);
cachedNodes = cell(1,0);
cachedWeights = cell(1,0);
pointRows = find(observations.kind == "point");
cellRows = find(observations.kind ~= "point");
if columns == 0
    status(:) = "noColumns";
elseif count > 0
    if ~isempty(cellRows) && 2*initialOrder > options.maxOrder
        error("IMVerticalSampling:IntegrationBudgetExceeded","maxOrder=%d cannot compare the initial order %d with twice that order. Increase the budget or explicitly select a justified initialOrder.",options.maxOrder,initialOrder);
    end
    pointKeys = variables(pointRows)+":"+string(orders(pointRows));
    for key = reshape(unique(pointKeys,"stable"),1,[])
        selected = pointRows(pointKeys == key);
        for first = 1:1024:numel(selected)
            batch = selected(first:min(first+1023,numel(selected)));
            synthesis(batch,:) = evaluatePoints(observations.zLower(batch),variables(batch(1)),orders(batch(1)));
        end
        evaluationPointCount(selected) = 1;
    end
    for row = reshape(cellRows,1,[])
        lower = observations.zLower(row);
        upper = observations.zUpper(row);
        pieces = [lower;breakpoints(breakpoints > lower & breakpoints < upper);upper];
        subintervalCount(row) = numel(pieces)-1;
        divisor = 1;
        if observations.kind(row) == "cellAverage", divisor = upper-lower; end
        order = initialOrder;
        while 2*order <= options.maxOrder
            value = zeros(1,columns);
            estimate = zeros(1,columns);
            scale = zeros(1,columns);
            seenNonzero = false(1,columns);
            for piece = 1:numel(pieces)-1
                a = pieces(piece); b = pieces(piece+1);
                [coarse,~] = integratePiece(a,b,order,variables(row),orders(row),divisor);
                [fine,magnitude,nonzero] = integratePiece(a,b,2*order,variables(row),orders(row),divisor);
                value = value+fine;
                estimate = estimate+abs(fine-coarse);
                scale = scale+magnitude;
                seenNonzero = seenNonzero | nonzero;
            end
            evaluationPointCount(row) = evaluationPointCount(row)+3*order*subintervalCount(row);
            if any(~isfinite([value estimate scale]))
                error("IMVerticalSampling:UnrepresentableIntegral","Observation %s has a nonfinite integral, error estimate or L1 scale. Choose representable profile units and geometry.",observations.observationLabel(row));
            end
            if any(scale == 0 & seenNonzero)
                error("IMVerticalSampling:UnrepresentableIntegral","Observation %s lost all positive quadrature mass for a nonzero evaluated profile. Change depth units or profile normalization; a zero integral and convergence estimate would be misleading.",observations.observationLabel(row));
            end
            if all(estimate <= absoluteTolerance(row)+options.relativeTolerance*scale)
                synthesis(row,:) = value;
                estimatedError(row,:) = estimate;
                l1Scale(row,:) = scale;
                finalOrder(row) = 2*order;
                converged(row) = true;
                status(row) = "estimatedToleranceSatisfied";
                break
            end
            order = 2*order;
        end
        if ~converged(row)
            error("IMVerticalSampling:IntegrationBudgetExceeded","Observation %s did not satisfy the estimated integration tolerance within maxOrder=%d. Supply known breakpoints or increase the justified quadrature budget; no matrix is returned.",observations.observationLabel(row),options.maxOrder);
        end
    end
end
overallStatus = "notRequired";
if count > 0 && columns == 0
    overallStatus = "noColumns";
elseif ~isempty(cellRows)
    overallStatus = "estimatedToleranceSatisfied";
end
integration = struct(relativeTolerance=options.relativeTolerance,absoluteTolerance=absoluteTolerance,initialOrder=initialOrder,initialOrderPolicy=initialOrderPolicy,maxOrder=options.maxOrder,breakpoints=breakpoints,basisBreakpoints=info.breakpoints,callerBreakpoints=options.breakpoints,estimatedError=estimatedError,l1Scale=l1Scale,finalOrder=finalOrder,subintervalCount=subintervalCount,evaluationPointCount=evaluationPointCount,totalEvaluationPoints=sum(evaluationPointCount),converged=converged,status=status,overallStatus=overallStatus,rule="positive Gauss-Legendre refinement n versus 2n, with per-piece absolute differences",estimateStatus="estimated, not a certified error bound",sourceAccuracy="unverified",sourceHints=info);
record = struct(schemaVersion=1,observations=observations,observationLabels=self.observationLabels,columnLabels=basis.columnLabels,zDomain=basis.zDomain,scientificSpecification=basis.scientificSpecification,constructionRecord=basis.constructionRecord,derivativeConvention=basis.capabilities.derivativeConvention,traceConvention=basis.capabilities.traceConvention,resolvedVariables=variables,resolvedDerivativeOrders=orders,integration=integration);

    function values = evaluatePoints(z,variable,derivativeOrder)
        values = basis.evaluate(z,variable,derivativeOrder);
        if ~isequal(size(values),[numel(z) columns]) || any(~isfinite(values),"all")
            error("IMVerticalSampling:InvalidEvaluation","The declared continuous evaluator must return finite sample-by-column values.");
        end
    end

    function [value,magnitude,seenNonzero] = integratePiece(lower,upper,order,variable,derivativeOrder,divisor)
        [nodes,weights] = quadratureRule(order);
        width = upper-lower;
        midpoint = lower+width/2;
        if ~isfinite(width) || midpoint <= lower || midpoint >= upper
            error("IMVerticalSampling:InvalidSupport","Integration subdivisions require a finite width and a representable interior midpoint.");
        end
        value = zeros(1,columns);
        magnitude = zeros(1,columns);
        seenNonzero = false(1,columns);
        for firstNode = 1:1024:order
            nodeIndices = firstNode:min(firstNode+1023,order);
            z = lower+width*((nodes(nodeIndices)+1)/2);
            if any(z <= lower | z >= upper)
                error("IMVerticalSampling:InvalidSupport","Gauss nodes must remain strictly inside each integration piece in floating-point arithmetic. Use representable depth units or wider cells.");
            end
            values = evaluatePoints(z,variable,derivativeOrder);
            % For means, normalize the measure before multiplying profiles:
            % a representable mean need not have a representable raw integral.
            physicalWeights = ((width/divisor)/2)*weights(nodeIndices);
            if any(physicalWeights <= 0 | ~isfinite(physicalWeights))
                error("IMVerticalSampling:UnrepresentableIntegral","A positive physical quadrature measure became zero or nonfinite. Change depth units or observation geometry.");
            end
            value = value+physicalWeights.'*values;
            magnitude = magnitude+physicalWeights.'*abs(values);
            seenNonzero = seenNonzero | any(values ~= 0,1);
        end
    end

    function [nodes,weights] = quadratureRule(order)
        index = find(cachedOrders == order,1);
        if isempty(index)
            [nodes,weights] = legpts(order);
            nodes = nodes(:); weights = weights(:);
            if numel(nodes) ~= order || numel(weights) ~= order || any(~isfinite([nodes;weights])) || any(weights <= 0)
                error("IMVerticalSampling:InvalidQuadratureRule","A positive finite Gauss-Legendre rule is required.");
            end
            cachedOrders(end+1) = order;
            cachedNodes{end+1} = nodes;
            cachedWeights{end+1} = weights;
        else
            nodes = cachedNodes{index};
            weights = cachedWeights{index};
        end
    end
end
