classdef IMDiscreteTransform
    % Store forward and inverse matrices for a scalar modal transform.
    %
    % For sampled values $$\mathbf{x}$$ and retained modal coefficients
    % $$\mathbf{a}$$, the transform matrices satisfy
    %
    % $$
    % \mathbf{a}=A_{\mathrm f}\mathbf{x},
    % \qquad
    % \widehat{\mathbf{x}}=A_{\mathrm i}\mathbf{a},
    % $$
    %
    % where `forwardMatrix` is $$A_{\mathrm f}$$ and `inverseMatrix` is
    % $$A_{\mathrm i}=\Phi$$, the sampled modal basis. The Galerkin forward
    % matrix is
    %
    % $$
    % A_{\mathrm f}=(\Phi^\mathsf{T}W\Phi)^{-1}\Phi^\mathsf{T}W.
    % $$
    %
    % The matrices may be rectangular. They obey
    %
    % $$
    % A_{\mathrm f}A_{\mathrm i}=I,
    % \qquad
    % A_{\mathrm i}A_{\mathrm f}=P_W,
    % $$
    %
    % where $$P_W$$ is the $$W$$-Galerkin projector onto the retained
    % sampled modal subspace. When $$W$$ is positive definite, this is the
    % $$W$$-orthogonal projector.
    %
    % Construct transforms from solved scalar modes with
    % `IMBasisSet.discreteTransform`.
    %
    % ```matlab
    % transform = basisSet.discreteTransform(z=z,nModes=8);
    % coefficients = transform.project(values);
    % valuesFit = transform.reconstruct(coefficients);
    % ```
    %
    % - Topic: Create discrete transforms
    % - Topic: Use transform matrices
    % - Topic: Apply discrete transforms
    % - Topic: Inspect samples and modes
    % - Topic: Assess transform quality
    % - Topic: Developer topics
    % - Declaration: classdef IMDiscreteTransform

    properties (SetAccess = private)
        % Forward transform matrix $$A_{\mathrm f}$$.
        %
        % Multiplying sampled profiles by `forwardMatrix` returns retained
        % modal coefficients. This is the matrix applied by `project`.
        %
        % - Topic: Use transform matrices
        % - nav_order: 10
        forwardMatrix

        % Inverse transform matrix $$A_{\mathrm i}=\Phi$$.
        %
        % Multiplying modal coefficients by `inverseMatrix` reconstructs
        % profiles on `z`. The matrix is the sampled modal basis and may be
        % rectangular. This is the matrix applied by `reconstruct`.
        %
        % - Topic: Use transform matrices
        % - nav_order: 20
        inverseMatrix

        % Physical sample points.
        %
        % - Topic: Inspect samples and modes
        % - nav_order: 10
        z

        % Quadrature increments aligned with `z`.
        %
        % Signed increments are retained so algebraic quadrature rules can
        % be inspected explicitly.
        %
        % - Topic: Inspect samples and modes
        % - nav_order: 20
        increments

        % Retained physical mode labels.
        %
        % - Topic: Inspect samples and modes
        % - nav_order: 30
        modeNumber

        % Basis-set normalization used to sample the modes.
        %
        % - Topic: Inspect samples and modes
        % - nav_order: 40
        normalization

        % Signed-norm-scaled relative Gram error.
        %
        % With
        % $$S=\operatorname{diag}(|\operatorname{diag}\Gamma_0|^{-1/2})$$,
        % this is
        % $$\|S(\Gamma-\Gamma_0)S\|_2$$. When $$\Gamma_0$$ is positive
        % definite, this is the worst-case relative Parseval error.
        %
        % - Topic: Assess transform quality
        % - nav_order: 10
        relativeGramError

        % Coefficient round-trip error $$\|A_{\mathrm f}A_{\mathrm i}-I\|_2$$.
        %
        % - Topic: Assess transform quality
        % - nav_order: 20
        roundTripError

        % Two-norm condition number of `inverseMatrix`.
        %
        % - Topic: Assess transform quality
        % - nav_order: 30
        inverseMatrixConditionNumber

        % Two-norm condition number of the sampled Gram matrix.
        %
        % - Topic: Assess transform quality
        % - nav_order: 40
        gramConditionNumber

        % Whether every diagonal entry of $$\Gamma_0$$ is positive.
        %
        % - Topic: Assess transform quality
        % - nav_order: 50
        targetGramIsPositiveDefinite

        % Whether at least one supplied quadrature increment is negative.
        %
        % - Topic: Assess transform quality
        % - nav_order: 60
        hasNegativeIncrements

        % Sample-space metric matrix $$W$$.
        %
        % - Topic: Developer topics — Inspect metric construction
        % - nav_order: 10
        % - Developer: true
        metricMatrix

        % Sampled Gram matrix $$\Gamma=A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}$$.
        %
        % - Topic: Developer topics — Inspect metric construction
        % - nav_order: 20
        % - Developer: true
        gramMatrix

        % Continuous diagonal Gram target $$\Gamma_0$$.
        %
        % - Topic: Developer topics — Inspect metric construction
        % - nav_order: 30
        % - Developer: true
        targetGramMatrix
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
            % - Parameter options.inverseMatrix: inverse transform matrix containing the sampled modes
            % - Parameter options.metricMatrix: sampled metric matrix
            % - Parameter options.targetGramMatrix: continuous diagonal Gram target
            % - Returns transform: initialized scalar discrete transform
            arguments
                options.z (:,1) double {mustBeReal, mustBeFinite}
                options.increments (:,1) double {mustBeReal, mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.normalization {mustBeTextScalar}
                options.inverseMatrix (:,:) double {mustBeReal, mustBeFinite}
                options.metricMatrix (:,:) double {mustBeReal, mustBeFinite}
                options.targetGramMatrix (:,:) double {mustBeReal, mustBeFinite}
            end

            z = options.z(:);
            increments = options.increments(:);
            inverseMatrix = options.inverseMatrix;
            metricMatrix = options.metricMatrix;
            targetGramMatrix = options.targetGramMatrix;
            nSamples = length(z);
            nModes = size(inverseMatrix,2);

            if length(increments) ~= nSamples || size(inverseMatrix,1) ~= nSamples
                error("IMDiscreteTransform:InvalidShape", "z, increments, and inverseMatrix rows must describe the same sample count.");
            end
            if length(options.modeNumber) ~= nModes
                error("IMDiscreteTransform:InvalidShape", "modeNumber must contain one label for each inverseMatrix column.");
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
            gramMatrix = inverseMatrix.'*metricMatrix*inverseMatrix;
            gramMatrix = 0.5*(gramMatrix + gramMatrix.');
            singularValues = svd(gramMatrix);
            rankTolerance = max(size(gramMatrix))*eps(norm(gramMatrix,2));
            if sum(singularValues > rankTolerance) < nModes
                error("IMDiscreteTransform:SingularGramMatrix", "The sampled Gram matrix is numerically rank deficient for the requested modes.");
            end

            forwardMatrix = gramMatrix \ (inverseMatrix.'*metricMatrix);
            scale = 1./sqrt(abs(targetNorms));
            scaledDifference = scale.*(gramMatrix - targetGramMatrix).*scale.';

            self.z = z;
            self.increments = increments;
            self.modeNumber = reshape(options.modeNumber,1,[]);
            self.normalization = string(options.normalization);
            self.inverseMatrix = inverseMatrix;
            self.metricMatrix = metricMatrix;
            self.forwardMatrix = forwardMatrix;
            self.gramMatrix = gramMatrix;
            self.targetGramMatrix = targetGramMatrix;
            self.relativeGramError = norm(scaledDifference,2);
            self.roundTripError = norm(forwardMatrix*inverseMatrix - eye(nModes),2);
            self.inverseMatrixConditionNumber = cond(inverseMatrix);
            self.gramConditionNumber = cond(gramMatrix);
            self.targetGramIsPositiveDefinite = all(targetNorms > 0);
            self.hasNegativeIncrements = any(increments < 0);
        end

        function coefficients = project(self, values)
            % Return retained modal coefficients for sampled profiles.
            %
            % This method applies `forwardMatrix`:
            %
            % $$
            % \mathbf{a}=A_{\mathrm f}\mathbf{x}.
            % $$
            %
            % - Topic: Apply discrete transforms
            % - nav_order: 10
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
            % Return sampled profiles reconstructed from retained coefficients.
            %
            % This method applies `inverseMatrix`:
            %
            % $$
            % \widehat{\mathbf{x}}=A_{\mathrm i}\mathbf{a}.
            % $$
            %
            % - Topic: Apply discrete transforms
            % - nav_order: 20
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
            values = self.inverseMatrix*coefficients;
        end
    end
end
