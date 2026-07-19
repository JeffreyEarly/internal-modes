function [transform, weightFit] = discreteTransform(self, options)
% Build a scalar Galerkin transform on fixed sample points.
%
% For retained normalized modes $$u_j$$ sampled at points $$z_i$$, this
% method forms
%
% $$
% (A_{\mathrm i})_{ij}=\Phi_{ij}=u_j(z_i),\qquad
% W_{\mathrm{int}}=\operatorname{diag}\!\left(r(z_i)w_i\right),
% $$
%
% then constructs the Galerkin forward matrix stored by
% `IMDiscreteTransform`. Eigenvalue-dependent endpoint terms are included
% when they depend only on a sampled endpoint value. Endpoint derivative
% traces cannot be inferred from arbitrary point samples and are rejected.
% When `weights` is omitted, `quadratureWeightsForPoints` chooses
% nonnegative weights with exact full-depth coverage by normalized Gram
% Frobenius fitting. This fitted path requires `lsqlin` from Optimization
% Toolbox.
%
% ```matlab
% [transform,weightFit] = basisSet.discreteTransform(z=z,nModes=8);
% transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
% coefficients = transform.transformForward(values);
% valuesFit = transform.transformBack(coefficients);
% ```
%
% When `weights` is omitted, the optional second output preserves the fit
% diagnostics and geometric comparison used to build the transform:
%
% ```matlab
% [weightFit.residualNorm weightFit.geometricResidualNorm]
% [weightFit.transform.relativeGramOperatorError weightFit.geometricTransform.relativeGramOperatorError]
% ```
%
% In this case `transform` is the same fitted transform stored in
% `weightFit.transform`. Supplying `weights` bypasses fitting, so the
% optional `weightFit` output is empty.
%
% - Topic: Build discrete transforms
% - Declaration: [transform,weightFit] = discreteTransform(basisSet,options)
% - Parameter options.z: increasing physical sample points
% - Parameter options.weights: optional quadrature weights aligned with `z`
% - Parameter options.nModes: number of leading retained modes
% - Returns transform: scalar discrete Galerkin transform
% - Returns weightFit: quadrature-fit diagnostics, or empty when weights are supplied
arguments
    self IMBasisSet
    options.z (:,1) double {mustBeReal, mustBeFinite}
    options.weights (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nModes (1,1) double {mustBeInteger, mustBePositive} = size(self.nativeModes,2)
end

z = options.z(:);
weights = options.weights(:);
nModes = options.nModes;

if isempty(weights)
    [~, weightFit] = self.quadratureWeightsForPoints(z=z,nModes=nModes);
    transform = weightFit.transform;
    return
end
weightFit = [];

if length(z) ~= length(weights)
    error("IMBasisSet:InvalidDiscreteWeights", "weights must contain one value for each sample point in z.");
end
if nModes > size(self.nativeModes,2)
    error("IMBasisSet:InvalidDiscreteModeCount", "The basis set contains %d modes, but nModes=%d was requested.", size(self.nativeModes,2), nModes);
end
if length(z) < nModes
    error("IMBasisSet:InsufficientDiscreteSamples", "At least %d sample points are required for %d retained modes.", nModes, nModes);
end
if length(z) < 2 || any(diff(z) <= 0)
    error("IMBasisSet:InvalidDiscreteGrid", "z must contain at least two strictly increasing, unique sample points.");
end
domainTolerance = 100*eps(max(1,max(abs(self.zDomain))));
if any(z < self.zDomain(1) - domainTolerance) || any(z > self.zDomain(2) + domainTolerance)
    error("IMBasisSet:InvalidDiscreteGrid", "All sample points must lie inside the basis-set zDomain.");
end
if ~any(weights ~= 0)
    error("IMBasisSet:InvalidDiscreteWeights", "weights must contain at least one nonzero value.");
end

inverseMatrix = self.u(z);
inverseMatrix = inverseMatrix(:,1:nModes);
context = self.evp.contextForSolver(self.solver);
spec = self.evp.innerProduct();
interiorWeight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight, z, context);
if isscalar(interiorWeight)
    interiorWeight = interiorWeight*ones(size(z));
elseif numel(interiorWeight) ~= length(z)
    error("IMBasisSet:InvalidDiscreteMetricWeight", "The EVP interior weight must return one value for each sample point.");
else
    interiorWeight = interiorWeight(:);
end
if ~isreal(interiorWeight) || any(~isfinite(interiorWeight))
    error("IMBasisSet:InvalidDiscreteMetricWeight", "The EVP interior weight must be finite and real at every sample point.");
end

metricMatrix = diag(interiorWeight.*weights);
if isfield(spec, "endpointInnerProductTerms") && ~isempty(spec.endpointInnerProductTerms)
    error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "Cross-variable endpoint inner products are not supported by scalar discrete transforms.");
end
endpointWeights = [spec.surfaceWeights; spec.bottomWeights];
for iWeight = 1:numel(endpointWeights)
    endpointWeight = endpointWeights(iWeight);
    if endpointWeight.d ~= 0
        error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "Endpoint derivative traces cannot be determined from arbitrary point samples; scalar discrete transforms require endpoint terms that depend only on sampled values.");
    end
    if endpointWeight.c == 0
        continue;
    end
    if string(endpointWeight.location) == "surface"
        zEndpoint = self.zDomain(2);
    else
        zEndpoint = self.zDomain(1);
    end
    endpointIndex = find(abs(z - zEndpoint) <= domainTolerance, 1);
    if isempty(endpointIndex)
        error("IMBasisSet:MissingDiscreteEndpointSample", "The %s endpoint must be included in z to represent its value-only inner-product term.", string(endpointWeight.location));
    end
    if ~isfinite(endpointWeight.coefficient)
        error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "The %s endpoint inner-product coefficient is not finite.", string(endpointWeight.location));
    end
    metricMatrix(endpointIndex,endpointIndex) = metricMatrix(endpointIndex,endpointIndex) ...
        + endpointWeight.coefficient*endpointWeight.c*endpointWeight.c;
end

continuousGramMatrix = self.gramMatrix();
targetGramMatrix = diag(diag(continuousGramMatrix(1:nModes,1:nModes)));
transform = IMDiscreteTransform(z=z,weights=weights,modeNumber=self.modeNumber(1:nModes),normalization=self.normalizationName(self.normalization),inverseMatrix=inverseMatrix,metricMatrix=metricMatrix,targetGramMatrix=targetGramMatrix);
end
