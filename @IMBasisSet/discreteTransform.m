function transform = discreteTransform(self, options)
% Build a scalar Galerkin transform on fixed sample points.
%
% For retained normalized modes $$u_j$$ sampled at points $$z_i$$, this
% method forms
%
% $$
% \Phi_{ij}=u_j(z_i),\qquad
% W_{\mathrm{int}}=\operatorname{diag}\!\left(r(z_i)\Delta z_i\right),
% $$
%
% then constructs the Galerkin forward matrix stored by
% `IMDiscreteTransform`. Eigenvalue-dependent endpoint terms are included
% when they depend only on a sampled endpoint value. Endpoint derivative
% traces cannot be inferred from arbitrary point samples and are rejected.
% When `increments` is omitted, `fitQuadrature` chooses nonnegative
% increments with exact full-depth coverage by normalized Gram fitting.
% This fitted path requires `lsqlin` from Optimization Toolbox.
%
% ```matlab
% transform = basisSet.discreteTransform(z=z,nModes=8);
% transform = basisSet.discreteTransform(z=z,increments=dz,nModes=8);
% coefficients = transform.project(values);
% ```
%
% - Topic: Build discrete transforms
% - Declaration: transform = discreteTransform(basisSet,options)
% - Parameter options.z: increasing physical sample points
% - Parameter options.increments: optional quadrature increments aligned with `z`
% - Parameter options.nModes: number of leading retained modes
% - Returns transform: scalar discrete Galerkin transform
arguments
    self IMBasisSet
    options.z (:,1) double {mustBeReal, mustBeFinite}
    options.increments (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nModes (1,1) double {mustBeInteger, mustBePositive} = size(self.nativeModes,2)
end

z = options.z(:);
increments = options.increments(:);
nModes = options.nModes;

if isempty(increments)
    fit = self.fitQuadrature(z=z, nModes=nModes);
    transform = fit.fittedTransform;
    return
end

if length(z) ~= length(increments)
    error("IMBasisSet:InvalidDiscreteIncrements", "increments must contain one value for each sample point in z.");
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
if ~any(increments ~= 0)
    error("IMBasisSet:InvalidDiscreteIncrements", "increments must contain at least one nonzero value.");
end

basisMatrix = self.u(z);
basisMatrix = basisMatrix(:,1:nModes);
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

metricMatrix = diag(interiorWeight.*increments);
if isfield(spec, "endpointInnerProductTerms") && ~isempty(spec.endpointInnerProductTerms)
    error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "Cross-variable endpoint inner products are not supported by scalar discrete transforms.");
end
endpointWeights = [spec.surfaceWeights; spec.bottomWeights];
for iWeight = 1:numel(endpointWeights)
    endpointWeight = endpointWeights(iWeight);
    if endpointWeight.d ~= 0
        error("IMBasisSet:UnsupportedDiscreteEndpointMetric", ...
            "Endpoint derivative traces cannot be determined from arbitrary point samples. Phase 1 supports value-only endpoint terms.");
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
        error("IMBasisSet:MissingDiscreteEndpointSample", ...
            "The %s endpoint must be included in z to represent its value-only inner-product term.", string(endpointWeight.location));
    end
    if ~isfinite(endpointWeight.coefficient)
        error("IMBasisSet:UnsupportedDiscreteEndpointMetric", "The %s endpoint inner-product coefficient is not finite.", string(endpointWeight.location));
    end
    metricMatrix(endpointIndex,endpointIndex) = metricMatrix(endpointIndex,endpointIndex) ...
        + endpointWeight.coefficient*endpointWeight.c*endpointWeight.c;
end

continuousGramMatrix = self.gramMatrix();
targetGramMatrix = diag(diag(continuousGramMatrix(1:nModes,1:nModes)));
transform = IMDiscreteTransform(z=z, increments=increments, modeNumber=self.modeNumber(1:nModes), ...
    normalization=self.normalizationName(self.normalization), basisMatrix=basisMatrix, ...
    metricMatrix=metricMatrix, targetGramMatrix=targetGramMatrix);
end
