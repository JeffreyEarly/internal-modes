function transform = scalarTransform(self,z,weights,options)
% Bind a genuine scalar physical projection to independently prepared recovery.
%
% The physical record retains the raw pairing P, continuous signed form H,
% positive majorant and normalized source dual B. Recovery uses the shared
% engine and applies R. The default correctedPhysicalDual policy corrects
% B only on the declared sampled span. A least-squares fit requires an
% explicit fitting or observation factor; quadrature is not a noise model.
%
% Fixed points and weights are preserved. This operation evaluates the
% existing continuous basis without another EVP solve or weight fit.
% Select or normalize the basis before preparing a new scalar transform.
%
% - Topic: Project modes
% - Parameter z: increasing physical sample column
% - Parameter weights: fixed physical dz quadrature weights
% - Parameter options.variable: supported scalar F/G channel, or APV for neutral labels
% - Parameter options.recovery: optional already prepared IMDiscreteTransform with matching synthesis, labels, activity and physical dual
% - Parameter options.policy: correctedPhysicalDual or sampledLeastSquares when preparing recovery
% - Parameter options.coefficientScale: nonsingular active coefficient scaling; the full continuous majorant defines the default
% - Parameter options.observationFactor: explicit positive observation whitening D
% - Parameter options.fitFactor: explicit possibly rectangular least-squares fitting factor
% - Parameter options.rankTolerance: dimensionless relative recovery rank cutoff
% - Parameter options.failureMode: diagnostic or error for unresolved recovery
% - Returns transform: IMScalarDiscreteTransform retaining independent physical and recovery records
arguments (Input)
    self (1,1) IMModeBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite}
    options.variable (1,1) string = ""
    options.recovery = []
    options.policy (1,1) string {mustBeMember(options.policy,["correctedPhysicalDual","sampledLeastSquares"])} = "correctedPhysicalDual"
    options.coefficientScale (:,:) double {mustBeFinite} = []
    options.observationFactor (:,:) double {mustBeFinite} = []
    options.fitFactor (:,:) double {mustBeFinite} = []
    options.rankTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-12
    options.failureMode (1,1) string {mustBeMember(options.failureMode,["diagnostic","error"])} = "diagnostic"
end
arguments (Output)
    transform (1,1) IMScalarDiscreteTransform
end
physical = self.projectionOnGrid(z,weights,variable=options.variable,prepareRecovery=false);
variable = options.variable;
if variable == ""
    if self.family == "neutralLabel", variable = "APV"; elseif self.family == "hydrostaticF", variable = "F"; else, variable = "G"; end
end
synthesis = self.evaluate(z,variable);
if isempty(options.recovery)
    if options.policy == "sampledLeastSquares" && isempty(options.fitFactor) && isempty(options.observationFactor)
        error("IMModeBasis:MissingFitFactor","sampledLeastSquares requires an explicit fitFactor or observationFactor; physical quadrature does not declare a fitting or noise norm.");
    end
    coefficientScale = options.coefficientScale;
    if isempty(coefficientScale)
        active = physical.activeColumnMask;
        majorant = physical.majorantGramMatrix(active,active);
        diagonalScale = 1./sqrt(diag(majorant));
        if isempty(diagonalScale)
            coefficientScale = zeros(0);
        else
            [factor,notPositive] = chol(diagonalScale.*majorant.*diagonalScale.');
            if notPositive ~= 0
                error("IMModeBasis:InvalidProjectionMajorant","The active continuous scalar majorant must be positive definite to define coefficient units.");
            end
            coefficientScale = diag(diagonalScale)/factor;
        end
    end
    recovery = IMDiscreteTransform.fromObservationMatrices(synthesis,policy=options.policy,columnLabels=self.columnLabels,activeColumnMask=physical.activeColumnMask,coefficientScale=coefficientScale,observationFactor=options.observationFactor,fitFactor=options.fitFactor,physicalDual=physical.physicalDualMatrix,rankTolerance=options.rankTolerance,failureMode=options.failureMode,metadata=struct(construction="IMModeBasis.scalarTransform",variable=variable,z=z,normalization=self.normalization,physicalProjectionProvenance=physical.provenance));
else
    if ~isa(options.recovery,"IMDiscreteTransform") || ~isscalar(options.recovery)
        error("IMModeBasis:InvalidPreparedRecovery","recovery must be one engine-prepared IMDiscreteTransform.");
    end
    if ~isempty(options.coefficientScale) || ~isempty(options.observationFactor) || ~isempty(options.fitFactor) || options.policy ~= "correctedPhysicalDual" || options.rankTolerance ~= 1e-12 || options.failureMode ~= "diagnostic"
        error("IMModeBasis:ConflictingRecoveryOptions","Supply either a prepared recovery or recovery preparation options.");
    end
    recovery = options.recovery;
end
if ~isequal(recovery.synthesisMatrix,synthesis)
    error("IMModeBasis:IncompatiblePreparedRecovery","The prepared recovery must use the exact requested scalar synthesis and observation order.");
end
if isfield(recovery.metadata,"variable") && string(recovery.metadata.variable) ~= variable
    error("IMModeBasis:IncompatiblePreparedRecovery","The prepared recovery variable must match the physical scalar channel.");
end
transform = IMScalarDiscreteTransform(physical,recovery,z=z,weights=weights,normalization=self.normalization);
end
