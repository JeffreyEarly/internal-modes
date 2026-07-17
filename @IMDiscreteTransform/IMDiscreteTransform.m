classdef IMDiscreteTransform
    % Store one scalar discrete Galerkin transform.
    %
    % `IMDiscreteTransform` stores a sampled basis $$\Phi$$, sampled metric
    % $$W$$, continuous diagonal target $$\Gamma_0$$, and the Galerkin
    % forward matrix
    %
    % $$
    % A_{\mathrm{gal}}=(\Phi^\mathsf{T}W\Phi)^{-1}\Phi^\mathsf{T}W.
    % $$
    %
    % Construct transforms from solved scalar modes with
    % `IMBasisSet.discreteTransform`.
    %
    % ```matlab
    % transform = basisSet.discreteTransform(z=z,increments=dz,nModes=8);
    % coefficients = transform.project(values);
    % valuesFit = transform.reconstruct(coefficients);
    % ```
    %
    % - Topic: Create discrete transforms
    % - Topic: Inspect discrete transforms
    % - Topic: Apply discrete transforms
    % - Declaration: classdef IMDiscreteTransform

    properties (SetAccess = private)
        % Physical sample points.
        %
        % - Topic: Inspect discrete transforms
        z

        % Quadrature increments aligned with `z`.
        %
        % Signed increments are retained so algebraic quadrature rules can
        % be inspected explicitly.
        %
        % - Topic: Inspect discrete transforms
        increments

        % Retained physical mode labels.
        %
        % - Topic: Inspect discrete transforms
        modeNumber

        % Basis-set normalization used to sample the modes.
        %
        % - Topic: Inspect discrete transforms
        normalization

        % Sampled reconstruction basis $$\Phi$$.
        %
        % Rows correspond to `z` and columns correspond to `modeNumber`.
        %
        % - Topic: Inspect discrete transforms
        basisMatrix

        % Sample-space metric matrix $$W$$.
        %
        % - Topic: Inspect discrete transforms
        metricMatrix

        % Galerkin forward matrix $$A_{\mathrm{gal}}$$.
        %
        % - Topic: Inspect discrete transforms
        forwardMatrix

        % Sampled Gram matrix $$\Gamma=\Phi^\mathsf{T}W\Phi$$.
        %
        % - Topic: Inspect discrete transforms
        gramMatrix

        % Continuous diagonal Gram target $$\Gamma_0$$.
        %
        % - Topic: Inspect discrete transforms
        targetGramMatrix

        % Signed-norm-scaled relative Gram error.
        %
        % With
        % $$S=\operatorname{diag}(|\operatorname{diag}\Gamma_0|^{-1/2})$$,
        % this is
        % $$\|S(\Gamma-\Gamma_0)S\|_2$$. When $$\Gamma_0$$ is positive
        % definite, this is the worst-case relative Parseval error.
        %
        % - Topic: Inspect discrete transforms
        relativeGramError

        % Coefficient round-trip error $$\|A_{\mathrm{gal}}\Phi-I\|_2$$.
        %
        % - Topic: Inspect discrete transforms
        roundTripError

        % Two-norm condition number of the sampled basis matrix.
        %
        % - Topic: Inspect discrete transforms
        basisConditionNumber

        % Two-norm condition number of the sampled Gram matrix.
        %
        % - Topic: Inspect discrete transforms
        gramConditionNumber

        % Whether every diagonal entry of $$\Gamma_0$$ is positive.
        %
        % - Topic: Inspect discrete transforms
        targetGramIsPositiveDefinite

        % Whether at least one supplied quadrature increment is negative.
        %
        % - Topic: Inspect discrete transforms
        hasNegativeIncrements
    end

    methods
        function self = IMDiscreteTransform(options)
            % Create a scalar discrete transform from canonical matrices.
            %
            % Ordinary users construct this object with
            % `IMBasisSet.discreteTransform`. This constructor is the
            % canonical matrix-level initialization path for later discrete
            % transform builders.
            %
            % - Topic: Create discrete transforms
            % - Declaration: transform = IMDiscreteTransform(options)
            % - Parameter options.z: physical sample points
            % - Parameter options.increments: quadrature increments
            % - Parameter options.modeNumber: retained mode labels
            % - Parameter options.normalization: basis normalization name
            % - Parameter options.basisMatrix: sampled basis matrix
            % - Parameter options.metricMatrix: sampled metric matrix
            % - Parameter options.targetGramMatrix: continuous diagonal Gram target
            % - Returns transform: initialized scalar discrete transform
            arguments
                options.z (:,1) double {mustBeReal, mustBeFinite}
                options.increments (:,1) double {mustBeReal, mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.normalization {mustBeTextScalar}
                options.basisMatrix (:,:) double {mustBeReal, mustBeFinite}
                options.metricMatrix (:,:) double {mustBeReal, mustBeFinite}
                options.targetGramMatrix (:,:) double {mustBeReal, mustBeFinite}
            end

            z = options.z(:);
            increments = options.increments(:);
            basisMatrix = options.basisMatrix;
            metricMatrix = options.metricMatrix;
            targetGramMatrix = options.targetGramMatrix;
            nSamples = length(z);
            nModes = size(basisMatrix,2);

            if length(increments) ~= nSamples || size(basisMatrix,1) ~= nSamples
                error("IMDiscreteTransform:InvalidShape", "z, increments, and basisMatrix rows must describe the same sample count.");
            end
            if length(options.modeNumber) ~= nModes
                error("IMDiscreteTransform:InvalidShape", "modeNumber must contain one label for each basisMatrix column.");
            end
            if ~isequal(size(metricMatrix), [nSamples nSamples])
                error("IMDiscreteTransform:InvalidShape", "metricMatrix must have one row and column for each sample point.");
            end
            if ~isequal(size(targetGramMatrix), [nModes nModes])
                error("IMDiscreteTransform:InvalidShape", "targetGramMatrix must have one row and column for each retained mode.");
            end

            metricTolerance = 100*eps(max(1,norm(metricMatrix,2)));
            if norm(metricMatrix - metricMatrix.',2) > metricTolerance
                error("IMDiscreteTransform:NonSymmetricMetric", "metricMatrix must be symmetric.");
            end
            targetTolerance = 100*eps(max(1,norm(targetGramMatrix,2)));
            if norm(targetGramMatrix - diag(diag(targetGramMatrix)),2) > targetTolerance
                error("IMDiscreteTransform:InvalidTargetGramMatrix", "targetGramMatrix must be diagonal.");
            end
            targetNorms = diag(targetGramMatrix);
            normTolerance = 100*eps(max(1,max(abs(targetNorms))));
            if any(~isfinite(targetNorms)) || any(abs(targetNorms) <= normTolerance)
                error("IMDiscreteTransform:InvalidTargetGramMatrix", "targetGramMatrix must have finite, numerically nonzero diagonal entries.");
            end

            metricMatrix = 0.5*(metricMatrix + metricMatrix.');
            gramMatrix = basisMatrix.'*metricMatrix*basisMatrix;
            gramMatrix = 0.5*(gramMatrix + gramMatrix.');
            singularValues = svd(gramMatrix);
            rankTolerance = max(size(gramMatrix))*eps(norm(gramMatrix,2));
            if sum(singularValues > rankTolerance) < nModes
                error("IMDiscreteTransform:SingularGramMatrix", "The sampled Gram matrix is numerically rank deficient for the requested modes.");
            end

            forwardMatrix = gramMatrix \ (basisMatrix.'*metricMatrix);
            scale = 1./sqrt(abs(targetNorms));
            scaledDifference = scale.*(gramMatrix - targetGramMatrix).*scale.';

            self.z = z;
            self.increments = increments;
            self.modeNumber = reshape(options.modeNumber,1,[]);
            self.normalization = string(options.normalization);
            self.basisMatrix = basisMatrix;
            self.metricMatrix = metricMatrix;
            self.forwardMatrix = forwardMatrix;
            self.gramMatrix = gramMatrix;
            self.targetGramMatrix = targetGramMatrix;
            self.relativeGramError = norm(scaledDifference,2);
            self.roundTripError = norm(forwardMatrix*basisMatrix - eye(nModes),2);
            self.basisConditionNumber = cond(basisMatrix);
            self.gramConditionNumber = cond(gramMatrix);
            self.targetGramIsPositiveDefinite = all(targetNorms > 0);
            self.hasNegativeIncrements = any(increments < 0);
        end

        function coefficients = project(self, values)
            % Project sampled profiles onto retained modal coefficients.
            %
            % - Topic: Apply discrete transforms
            % - Declaration: coefficients = project(transform,values)
            % - Parameter values: sampled profiles with rows aligned to `z`
            % - Returns coefficients: retained modal coefficients
            arguments
                self IMDiscreteTransform
                values double {mustBeFinite}
            end

            if size(values,1) ~= length(self.z)
                error("IMDiscreteTransform:InvalidSampleCount", "values must have one row for each sample point in z.");
            end
            coefficients = self.forwardMatrix*values;
        end

        function values = reconstruct(self, coefficients)
            % Reconstruct sampled profiles from retained coefficients.
            %
            % - Topic: Apply discrete transforms
            % - Declaration: values = reconstruct(transform,coefficients)
            % - Parameter coefficients: coefficient arrays with rows aligned to `modeNumber`
            % - Returns values: reconstructed profiles sampled on `z`
            arguments
                self IMDiscreteTransform
                coefficients double {mustBeFinite}
            end

            if size(coefficients,1) ~= length(self.modeNumber)
                error("IMDiscreteTransform:InvalidCoefficientCount", "coefficients must have one row for each retained mode.");
            end
            values = self.basisMatrix*coefficients;
        end
    end
end
