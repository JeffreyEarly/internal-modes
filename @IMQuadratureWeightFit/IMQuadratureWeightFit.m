classdef IMQuadratureWeightFit
    % Store diagnostics for quadrature weights fitted on fixed points.
    %
    % A quadrature-weight fit holds the physical sample points $$z_k$$, the
    % retained modes, and their normalization fixed. Only one quadrature
    % weight $$w_k$$ per point is optimized. For the sampled mode matrix
    %
    % $$
    % \Phi_{ki}=u_i(z_k),
    % $$
    %
    % the weights produce the discrete Gram matrix
    %
    % $$
    % \Gamma(w)
    % =\Phi^{\mathsf T}
    % \operatorname{diag}\!\left(r(z_k)w_k\right)\Phi
    % +\Gamma_{\mathrm{endpoint}}.
    % $$
    %
    % Here $$r(z)$$ is the EVP interior inner-product weight and
    % $$\Gamma_{\mathrm{endpoint}}$$ contains endpoint contributions that do
    % not depend on the interior quadrature weights. The continuous Gram
    % target $$\Gamma_0$$ contains the desired inner products of the retained
    % normalized modes. The default objective adjusts $$w$$ so that
    % $$\Gamma(w)\approx\Gamma_0$$. With
    % $$C_i=(\Gamma_0)_{ii}$$, its least-squares system is
    %
    % $$
    % (A_{\mathrm{LS}})_{(i,j),k}
    % =\frac{r(z_k)\Phi_{ki}\Phi_{kj}}{\sqrt{|C_iC_j|}},
    % \qquad
    % (b_{\mathrm{LS}})_{(i,j)}
    % =\frac{(\Gamma_0-\Gamma_{\mathrm{endpoint}})_{ij}}
    % {\sqrt{|C_iC_j|}},
    % $$
    %
    % and the fitted weights solve
    %
    % $$
    % \min_w\left\|A_{\mathrm{LS}}w-b_{\mathrm{LS}}\right\|_2.
    % $$
    %
    % By default the fit also requires
    %
    % $$
    % w_k\geq0,
    % \qquad
    % \sum_k w_k=D,
    % $$
    %
    % where $$D=z_s-z_b$$ is the full physical depth.
    %
    % Before optimization, the fixed points define a natural geometric
    % control-volume rule. For $$n$$ points, let
    %
    % $$
    % e_1=z_b,
    % \qquad
    % e_{k+1}=\frac{z_k+z_{k+1}}{2}\quad(k=1,\ldots,n-1),
    % \qquad
    % e_{n+1}=z_s.
    % $$
    %
    % The unoptimized weights are
    %
    % $$
    % w_k^{\mathrm{geometric}}=e_{k+1}-e_k.
    % $$
    %
    % These are the literal midpoint/control-volume approximation to the
    % sampled integral: they use only point geometry and do not use modal
    % orthogonality. `quadratureWeightsForPoints` supplies them as the
    % optimizer's initial guess. They are positive, cover the full domain,
    % and define `geometricTransform`, the reference obtained by skipping
    % optimization. The fitted algebraic weights are instead adjusted to
    % reproduce the retained modal inner products and need not retain a
    % geometric-width interpretation when constraints are relaxed.
    %
    % `transform` and `weights` are the optimized production result;
    % `geometricTransform` and `geometricWeights` are the unoptimized
    % comparison. `residualNorm` and `geometricResidualNorm` evaluate the
    % same fitting objective for those two rules. Their transforms'
    % `relativeGramError` values measure the resulting Gram mismatch.
    % Constraint properties record whether nonnegativity and full-depth
    % coverage were imposed, while `depthError` and
    % `transform.hasNegativeWeights` report the corresponding fitted result.
    % A custom objective changes $$A_{\mathrm{LS}}$$ and
    % $$b_{\mathrm{LS}}$$, but the fitted and geometric rules still use the
    % same points, modes, and normalization.
    %
    % Obtain the fitted weights as the primary output and this diagnostic
    % object as the optional second output of
    % `IMBasisSet.quadratureWeightsForPoints`. The returned `weights`,
    % `weightFit.weights`, and `weightFit.transform.weights` are the same
    % vector.
    %
    % ```matlab
    % [weights,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
    % weightFit.residualNorm
    % weightFit.geometricResidualNorm
    % weightFit.transform.relativeGramError
    % weightFit.geometricTransform.relativeGramError
    % coefficients = weightFit.transform.transformForward(values);
    % ```
    %
    % - Topic: Inspect fitted weights
    % - Topic: Compare geometric weights
    % - Topic: Assess fit quality
    % - Topic: Inspect constraints
    % - Topic: Developer topics
    % - Declaration: classdef IMQuadratureWeightFit

    properties (SetAccess = private)
        % Discrete transform constructed with the fitted weights.
        %
        % This is the production transform associated with `weights`. Its
        % sample points, retained modes, normalization, metric, and transform
        % matrices match the quadrature-weight fit.
        %
        % ```matlab
        % coefficients = weightFit.transform.transformForward(values);
        % ```
        %
        % - Topic: Inspect fitted weights
        % - nav_order: 10
        transform

        % Fitted algebraic quadrature weights aligned with the fixed points.
        %
        % These are the weights $$w_k$$ returned as the first output of
        % `quadratureWeightsForPoints`. They minimize the selected objective
        % subject to the requested constraints and satisfy
        % `weights = transform.weights`.
        %
        % - Topic: Inspect fitted weights
        % - nav_order: 20
        weights

        % Reference transform using geometric control-volume weights.
        %
        % This transform uses the same points, modes, and normalization as
        % `transform`, but its weights come only from the physical control
        % volumes around the fixed points. It provides a baseline for Gram
        % and transform-quality comparisons.
        %
        % - Topic: Compare geometric weights
        % - nav_order: 10
        geometricTransform

        % Geometric control-volume weights aligned with the fixed points.
        %
        % For adjacent points, control-volume edges lie at their midpoints;
        % the outer edges are the basis-set domain boundaries. Therefore
        % these weights are geometric widths and sum to the full depth.
        % They satisfy `geometricWeights = geometricTransform.weights`.
        %
        % - Topic: Compare geometric weights
        % - nav_order: 20
        geometricWeights

        % Name of the least-squares objective used to fit the weights.
        %
        % The default value is `"normalizedGram"`. A custom objective may
        % provide its own name in the returned specification struct.
        %
        % - Topic: Assess fit quality
        % - nav_order: 10
        objectiveName

        % Two-norm of the fitted objective residual.
        %
        % This is
        %
        % $$
        % \left\|A_{\mathrm{LS}}w-b_{\mathrm{LS}}\right\|_2.
        % $$
        %
        % Compare it with `geometricResidualNorm` to assess improvement over
        % the geometric control-volume weights under the same objective.
        %
        % - Topic: Assess fit quality
        % - nav_order: 20
        residualNorm

        % Two-norm of the geometric-weight objective residual.
        %
        % This is
        %
        % $$
        % \left\|A_{\mathrm{LS}}w_{\mathrm{geometric}}
        % -b_{\mathrm{LS}}\right\|_2.
        % $$
        %
        % - Topic: Assess fit quality
        % - nav_order: 30
        geometricResidualNorm

        % Whether the optimization required nonnegative weights.
        %
        % When true, the least-squares problem imposed $$w_k\geq0$$ at every
        % fixed point. When false, `transform.hasNegativeWeights` reports
        % whether the fitted algebraic rule actually contains negative
        % weights.
        %
        % - Topic: Inspect constraints
        % - nav_order: 10
        nonnegativeConstraint

        % Whether the weights were constrained to cover the full depth.
        %
        % When true, the least-squares problem imposed
        % $$\sum_k w_k=D$$, where $$D$$ is `depthTarget`.
        %
        % - Topic: Inspect constraints
        % - nav_order: 20
        depthConstraint

        % Full physical depth targeted by the weight-sum constraint.
        %
        % This is
        % $$D=z_\mathrm{surface}-z_\mathrm{bottom}$$ in the same coordinate
        % units as the fixed points, normally meters.
        %
        % - Topic: Inspect constraints
        % - nav_order: 30
        depthTarget

        % Difference between the fitted weight sum and `depthTarget`.
        %
        % The value is
        %
        % $$
        % \sum_k w_k-D.
        % $$
        %
        % It should be near zero when `depthConstraint` is true.
        %
        % - Topic: Inspect constraints
        % - nav_order: 40
        depthError

        % Difference between the geometric weight sum and `depthTarget`.
        %
        % The geometric control-volume construction covers the full domain,
        % so this value should be near roundoff.
        %
        % - Topic: Inspect constraints
        % - nav_order: 50
        geometricDepthError

        % Least-squares matrix $$A_{\mathrm{LS}}$$.
        %
        % It has one column per fixed point. Multiplying this matrix by a
        % physical quadrature-weight vector produces the modeled objective
        % quantities before subtracting `objectiveTarget`.
        %
        % - Topic: Developer topics — Inspect least-squares system
        % - nav_order: 10
        % - Developer: true
        objectiveMatrix

        % Least-squares target vector $$b_{\mathrm{LS}}$$.
        %
        % Together with `objectiveMatrix`, this defines the physical-unit
        % residual $$A_{\mathrm{LS}}w-b_{\mathrm{LS}}$$.
        %
        % - Topic: Developer topics — Inspect least-squares system
        % - nav_order: 20
        % - Developer: true
        objectiveTarget

        % Objective residual vector at the fitted weights.
        %
        % This is `objectiveMatrix*weights-objectiveTarget`.
        %
        % - Topic: Developer topics — Inspect least-squares system
        % - nav_order: 30
        % - Developer: true
        residual

        % Objective residual vector at the geometric weights.
        %
        % This is `objectiveMatrix*geometricWeights-objectiveTarget`.
        %
        % - Topic: Developer topics — Inspect least-squares system
        % - nav_order: 40
        % - Developer: true
        geometricResidual

        % Exit flag returned by `lsqlin`.
        %
        % Positive values indicate successful termination. Consult MATLAB's
        % `lsqlin` documentation together with `solverOutput` for detailed
        % optimizer diagnostics.
        %
        % - Topic: Developer topics — Inspect optimizer output
        % - nav_order: 10
        % - Developer: true
        exitFlag

        % Solver diagnostics returned by `lsqlin`.
        %
        % The struct contains the optimizer iteration count, termination
        % message, and other release-dependent `lsqlin` diagnostics.
        %
        % - Topic: Developer topics — Inspect optimizer output
        % - nav_order: 20
        % - Developer: true
        solverOutput
    end

    methods
        function self = IMQuadratureWeightFit(options)
            % Create diagnostics for a quadrature-weight fit.
            %
            % This developer constructor combines transforms built on the
            % same fixed points and retained modes with the least-squares
            % system and optimizer result that produced the fitted weights.
            % Ordinary users obtain this object from
            % `IMBasisSet.quadratureWeightsForPoints`.
            %
            % - Topic: Developer topics — Construct fit results
            % - Declaration: weightFit = IMQuadratureWeightFit(options)
            % - Parameter options.transform: transform using fitted weights
            % - Parameter options.geometricTransform: reference transform using geometric weights
            % - Parameter options.objectiveName: least-squares objective name
            % - Parameter options.objectiveMatrix: least-squares matrix
            % - Parameter options.objectiveTarget: least-squares target
            % - Parameter options.nonnegativeConstraint: whether nonnegative weights were imposed
            % - Parameter options.depthConstraint: whether exact depth was imposed
            % - Parameter options.depthTarget: full physical depth
            % - Parameter options.exitFlag: optimizer exit flag
            % - Parameter options.solverOutput: optimizer diagnostics
            % - Returns weightFit: initialized quadrature-weight fit diagnostics
            % - Developer: true
            arguments
                options.transform (1,1) IMDiscreteTransform
                options.geometricTransform (1,1) IMDiscreteTransform
                options.objectiveName {mustBeTextScalar}
                options.objectiveMatrix (:,:) double {mustBeReal, mustBeFinite}
                options.objectiveTarget (:,1) double {mustBeReal, mustBeFinite}
                options.nonnegativeConstraint (1,1) logical
                options.depthConstraint (1,1) logical
                options.depthTarget (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.exitFlag (1,1) double {mustBeReal, mustBeFinite}
                options.solverOutput (1,1) struct
            end

            transform = options.transform;
            geometricTransform = options.geometricTransform;
            if ~isequal(transform.z,geometricTransform.z) || ~isequal(transform.modeNumber,geometricTransform.modeNumber) || transform.normalization ~= geometricTransform.normalization
                error("IMQuadratureWeightFit:IncompatibleTransforms", "The fitted and geometric transforms must use the same points, modes, and normalization.");
            end
            if size(options.objectiveMatrix,2) ~= length(transform.z) || size(options.objectiveMatrix,1) ~= length(options.objectiveTarget)
                error("IMQuadratureWeightFit:InvalidObjective", "The objective dimensions must match the transform sample count and target length.");
            end

            residual = options.objectiveMatrix*transform.weights - options.objectiveTarget;
            geometricResidual = options.objectiveMatrix*geometricTransform.weights - options.objectiveTarget;

            self.transform = transform;
            self.weights = transform.weights;
            self.geometricTransform = geometricTransform;
            self.geometricWeights = geometricTransform.weights;
            self.objectiveName = string(options.objectiveName);
            self.objectiveMatrix = options.objectiveMatrix;
            self.objectiveTarget = options.objectiveTarget;
            self.residual = residual;
            self.geometricResidual = geometricResidual;
            self.residualNorm = norm(residual,2);
            self.geometricResidualNorm = norm(geometricResidual,2);
            self.nonnegativeConstraint = options.nonnegativeConstraint;
            self.depthConstraint = options.depthConstraint;
            self.depthTarget = options.depthTarget;
            self.depthError = sum(transform.weights) - options.depthTarget;
            self.geometricDepthError = sum(geometricTransform.weights) - options.depthTarget;
            self.exitFlag = options.exitFlag;
            self.solverOutput = options.solverOutput;
        end
    end
end
