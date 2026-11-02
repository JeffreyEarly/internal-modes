classdef IMScalarDiscreteTransform < IMDiscreteTransform
    % Attach scalar-grid geometry and physical projection data to recovery.
    %
    % The generic superclass owns synthesis and recovery. This adapter adds
    % scalar samples, quadrature and the genuine physical IMProjection record.
    % Its sampled Gram and continuous target remain distinct from numerical
    % sample rank. Failed diagnostic preparation retains synthesis/physical
    % pairings but does not expose an applicable recovery matrix.
    %
    % ```matlab
    % transform = basisSet.discreteTransform(z=z,weights=w,nModes=8);
    % physicalCoefficients = transform.projection.projectPhysical(values);
    % recoveredCoefficients = transform.transformForward(values);
    % ```
    %
    % - Topic: Create scalar transforms
    % - Topic: Inspect scalar transforms
    % - Declaration: classdef IMScalarDiscreteTransform < IMDiscreteTransform

    properties (SetAccess = private)
        % Physical sample points, one per scalar observation.
        % - Topic: Inspect scalar transforms
        z (:,1) double
        % Physical quadrature weights, distinct from observation-noise weights.
        % - Topic: Inspect scalar transforms
        weights (:,1) double
        % Scientific scalar mode numbers aligned with columnLabels.
        % - Topic: Inspect scalar transforms
        modeNumber (1,:) double
        % Continuous normalization captured by the scalar family.
        % - Topic: Inspect scalar transforms
        normalization (1,1) string
        % Physical pairing, continuous forms and independent Gram diagnostics.
        % - Topic: Inspect scalar transforms
        projection (1,1) IMProjection
    end

    properties (Dependent, SetAccess = private)
        % Sampled scalar synthesis, also available as synthesisMatrix.
        % - Topic: Inspect scalar transforms
        inverseMatrix
        % Genuine scalar physical sample metric.
        % - Topic: Inspect scalar transforms
        metricMatrix
        % Sampled physical Gram, not the recovery overlap or sample rank.
        % - Topic: Inspect scalar transforms
        gramMatrix
        % Continuous scalar physical Gram target.
        % - Topic: Inspect scalar transforms
        targetGramMatrix
        % Magnitude-normalized physical Gram discrepancy.
        % - Topic: Inspect scalar transforms
        relativeGramOperatorError
        % Numerical recovery round trip, infinity when recovery failed.
        % - Topic: Inspect scalar transforms
        roundTripError
        % Unscaled scalar synthesis condition number.
        % - Topic: Inspect scalar transforms
        inverseMatrixConditionNumber
        % Condition number of the actual sampled physical Gram.
        % - Topic: Inspect scalar transforms
        gramConditionNumber
        % Rank of the actual sampled physical Gram.
        % - Topic: Inspect scalar transforms
        sampledGramRank
        % Whether the continuous active scalar target is positive definite.
        % - Topic: Inspect scalar transforms
        targetGramIsPositiveDefinite
        % Whether any scalar physical quadrature weight is negative.
        % - Topic: Inspect scalar transforms
        hasNegativeWeights
    end

    methods
        function self = IMScalarDiscreteTransform(options)
            % Prepare one scalar adapter from its genuine physical matrices.
            %
            % IMProjection delegates recovery to the shared kernel. Its
            % prepared generic transform initializes the superclass without
            % another factorization. The physical dual remains separately
            % available through projection even when recovery fails.
            %
            % - Topic: Create scalar transforms
            % - Parameter options.z: scalar sample column
            % - Parameter options.weights: physical quadrature weights
            % - Parameter options.modeNumber: aligned scientific scalar labels
            % - Parameter options.normalization: continuous normalization name
            % - Parameter options.inverseMatrix: scalar sampled synthesis
            % - Parameter options.metricMatrix: symmetric physical sample metric
            % - Parameter options.targetGramMatrix: diagonal continuous scalar Gram
            % - Parameter options.failureMode: diagnostic or error on unresolved recovery
            % - Returns self: scalar adapter with generic recovery capability
            arguments (Input)
                options.z (:,1) double {mustBeReal,mustBeFinite}
                options.weights (:,1) double {mustBeReal,mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.normalization (1,1) string
                options.inverseMatrix (:,:) double {mustBeReal,mustBeFinite}
                options.metricMatrix (:,:) double {mustBeReal,mustBeFinite}
                options.targetGramMatrix (:,:) double {mustBeReal,mustBeFinite}
                options.failureMode (1,1) string = "diagnostic"
            end
            nSamples = numel(options.z); nModes = numel(options.modeNumber);
            if numel(options.weights) ~= nSamples || ~isequal(size(options.inverseMatrix),[nSamples nModes])
                error("IMScalarDiscreteTransform:InvalidShape","Samples, weights and mode labels must align with scalar synthesis.");
            end
            if ~isequal(size(options.targetGramMatrix),[nModes nModes]) || ~isequal(options.targetGramMatrix,diag(diag(options.targetGramMatrix)))
                error("IMScalarDiscreteTransform:InvalidTargetGram","Scalar targets must be diagonal and aligned with modeNumber.");
            end
            physical = IMProjection(options.inverseMatrix,options.metricMatrix,options.targetGramMatrix,columnLabels=string(options.modeNumber),failureMode=options.failureMode);
            self@IMDiscreteTransform(physical.recovery);
            self.z = options.z;
            self.weights = options.weights;
            self.modeNumber = options.modeNumber;
            self.normalization = options.normalization;
            self.projection = physical;
        end

        function value = get.inverseMatrix(self), value = self.synthesisMatrix; end
        function value = get.metricMatrix(self), value = self.projection.metricMatrix; end
        function value = get.gramMatrix(self), value = self.projection.gramMatrix; end
        function value = get.targetGramMatrix(self), value = self.projection.targetGramMatrix; end
        function value = get.relativeGramOperatorError(self), value = self.projection.gramError; end
        function value = get.roundTripError(self), value = self.projection.roundTripError; end
        function value = get.inverseMatrixConditionNumber(self), value = self.projection.inverseMatrixConditionNumber; end
        function value = get.gramConditionNumber(self), value = self.projection.gramConditionNumber; end
        function value = get.sampledGramRank(self), value = self.projection.sampledGramRank; end
        function value = get.targetGramIsPositiveDefinite(self), value = self.projection.targetGramIsPositiveDefinite; end
        function value = get.hasNegativeWeights(self), value = any(self.weights < 0); end
    end
end
