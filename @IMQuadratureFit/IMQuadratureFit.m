classdef IMQuadratureFit
    % Store a fixed-point quadrature fit and geometric comparison.
    %
    % `IMQuadratureFit` compares two scalar Galerkin transforms on the same
    % physical points: one using fitted increments and one using geometric
    % control-volume increments. If the least-squares objective is
    % $$\|A\Delta z-b\|_2$$, the stored residuals are
    % $$A\Delta z_\mathrm{fit}-b$$ and
    % $$A\Delta z_\mathrm{geometric}-b$$.
    %
    % Construct fits with `IMBasisSet.fitQuadrature`.
    %
    % ```matlab
    % fit = basisSet.fitQuadrature(z=z,nModes=8);
    % transform = fit.fittedTransform;
    % ```
    %
    % - Topic: Create quadrature fits
    % - Topic: Inspect quadrature fits
    % - Declaration: classdef IMQuadratureFit

    properties (SetAccess = private)
        % Transform using the fitted quadrature increments.
        %
        % - Topic: Inspect quadrature fits
        fittedTransform

        % Transform using geometric control-volume increments.
        %
        % - Topic: Inspect quadrature fits
        geometricTransform

        % Fitted quadrature increments.
        %
        % - Topic: Inspect quadrature fits
        fittedIncrements

        % Geometric control-volume increments.
        %
        % - Topic: Inspect quadrature fits
        geometricIncrements

        % Name of the least-squares objective.
        %
        % - Topic: Inspect quadrature fits
        objectiveName

        % Least-squares matrix $$A$$.
        %
        % - Topic: Inspect quadrature fits
        objectiveMatrix

        % Least-squares target $$b$$.
        %
        % - Topic: Inspect quadrature fits
        objectiveTarget

        % Objective residual at the fitted increments.
        %
        % - Topic: Inspect quadrature fits
        fittedObjectiveResidual

        % Objective residual at the geometric increments.
        %
        % - Topic: Inspect quadrature fits
        geometricObjectiveResidual

        % Two-norm of `fittedObjectiveResidual`.
        %
        % - Topic: Inspect quadrature fits
        fittedResidualNorm

        % Two-norm of `geometricObjectiveResidual`.
        %
        % - Topic: Inspect quadrature fits
        geometricResidualNorm

        % Whether the fit imposed nonnegative increments.
        %
        % - Topic: Inspect quadrature fits
        nonnegativeConstraint

        % Whether the fit imposed exact full-depth coverage.
        %
        % - Topic: Inspect quadrature fits
        depthConstraint

        % Requested full physical depth.
        %
        % - Topic: Inspect quadrature fits
        depthTarget

        % Difference between fitted increment sum and `depthTarget`.
        %
        % - Topic: Inspect quadrature fits
        fittedDepthError

        % Difference between geometric increment sum and `depthTarget`.
        %
        % - Topic: Inspect quadrature fits
        geometricDepthError

        % Exit flag returned by `lsqlin`.
        %
        % - Topic: Inspect quadrature fits
        exitFlag

        % Solver diagnostics returned by `lsqlin`.
        %
        % - Topic: Inspect quadrature fits
        solverOutput
    end

    methods
        function self = IMQuadratureFit(options)
            % Create a fixed-point quadrature fit result.
            %
            % Ordinary users construct this object with
            % `IMBasisSet.fitQuadrature`.
            %
            % - Topic: Create quadrature fits
            % - Declaration: fit = IMQuadratureFit(options)
            % - Parameter options.fittedTransform: transform using fitted increments
            % - Parameter options.geometricTransform: transform using geometric increments
            % - Parameter options.objectiveName: least-squares objective name
            % - Parameter options.objectiveMatrix: least-squares matrix
            % - Parameter options.objectiveTarget: least-squares target
            % - Parameter options.nonnegativeConstraint: whether nonnegative increments were imposed
            % - Parameter options.depthConstraint: whether exact depth was imposed
            % - Parameter options.depthTarget: full physical depth
            % - Parameter options.exitFlag: optimizer exit flag
            % - Parameter options.solverOutput: optimizer diagnostics
            % - Returns fit: initialized quadrature fit result
            arguments
                options.fittedTransform (1,1) IMDiscreteTransform
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

            fittedTransform = options.fittedTransform;
            geometricTransform = options.geometricTransform;
            if ~isequal(fittedTransform.z,geometricTransform.z) ...
                    || ~isequal(fittedTransform.modeNumber,geometricTransform.modeNumber) ...
                    || fittedTransform.normalization ~= geometricTransform.normalization
                error("IMQuadratureFit:IncompatibleTransforms", "The fitted and geometric transforms must use the same points, modes, and normalization.");
            end
            if size(options.objectiveMatrix,2) ~= length(fittedTransform.z) ...
                    || size(options.objectiveMatrix,1) ~= length(options.objectiveTarget)
                error("IMQuadratureFit:InvalidObjective", "The objective dimensions must match the transform sample count and target length.");
            end

            fittedResidual = options.objectiveMatrix*fittedTransform.increments - options.objectiveTarget;
            geometricResidual = options.objectiveMatrix*geometricTransform.increments - options.objectiveTarget;

            self.fittedTransform = fittedTransform;
            self.geometricTransform = geometricTransform;
            self.fittedIncrements = fittedTransform.increments;
            self.geometricIncrements = geometricTransform.increments;
            self.objectiveName = string(options.objectiveName);
            self.objectiveMatrix = options.objectiveMatrix;
            self.objectiveTarget = options.objectiveTarget;
            self.fittedObjectiveResidual = fittedResidual;
            self.geometricObjectiveResidual = geometricResidual;
            self.fittedResidualNorm = norm(fittedResidual,2);
            self.geometricResidualNorm = norm(geometricResidual,2);
            self.nonnegativeConstraint = options.nonnegativeConstraint;
            self.depthConstraint = options.depthConstraint;
            self.depthTarget = options.depthTarget;
            self.fittedDepthError = sum(fittedTransform.increments) - options.depthTarget;
            self.geometricDepthError = sum(geometricTransform.increments) - options.depthTarget;
            self.exitFlag = options.exitFlag;
            self.solverOutput = options.solverOutput;
        end
    end
end
