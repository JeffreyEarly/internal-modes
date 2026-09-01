function [transform, assessment] = certifiedDiscreteTransform(self, options)
% Select and independently fit a certified aligned F/G family band.
%
% Fresh family-specific weights are fitted for every count considered by
% the Gram search. This prevents a poorly resolved large candidate family
% from contaminating the assessment of a smaller MDA or APV band. Optional
% leakage and coupled quadratic policies can reduce that Gram-certified
% band; each reduction is refitted until the fitted and retained counts
% agree.
%
% - Topic: Build discrete transforms
% - Declaration: [transform,assessment] = certifiedDiscreteTransform(basisSet,options)
% - Parameter options.nPoints: exact requested mode-root point count
% - Parameter options.z: explicit increasing physical sample points
% - Parameter options.variables: requested direct channels, F and/or G
% - Parameter options.gridDesign: optional provenance returned by `modeRootGrid`
% - Parameter options.gramTolerance: per-channel normalized-Gram tolerance
% - Parameter options.leakageTolerance: optional rejected-mode leakage tolerance
% - Parameter options.quadraticAliasingTolerance: optional coupled-product tolerance
% - Parameter options.nCheckModes: optional rejected-mode check count
% - Returns transform: independently fitted certified aligned transform
% - Returns assessment: final diagnostics, grid provenance, and count-search table
arguments
    self IMInternalModesBasis
    options.nPoints double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.z (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.variables (1,:) string = strings(1,0)
    options.gridDesign struct = struct.empty
    options.gramTolerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-2
    options.leakageTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.quadraticAliasingTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nCheckModes double {mustBeReal, mustBeFinite} = zeros(0,1)
end

IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nPoints,"nPoints");
IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nCheckModes,"nCheckModes");
IMDiscreteTransformTools.validateOptionalPositiveTolerance(options.leakageTolerance,"leakageTolerance");
IMDiscreteTransformTools.validateOptionalPositiveTolerance(options.quadraticAliasingTolerance,"quadraticAliasingTolerance");
if self.evp.modeFamily == "meanDensityAnomaly" && ~isempty(options.quadraticAliasingTolerance)
    error("IMMeanDensityAnomalyModesBasis:UnavailableQuadraticAliasingPolicy", ...
        "Coupled quadratic aliasing is not defined for mean-density-anomaly modes.");
end
hasPointCount = ~isempty(options.nPoints);
hasExplicitPoints = ~isempty(options.z);
if hasPointCount == hasExplicitPoints
    error("IMBasisSet:InvalidDiscretePointSpecification", "Specify exactly one of nPoints or z.");
end
if hasPointCount && ~isempty(options.gridDesign)
    error("IMBasisSet:InvalidDiscretePointSpecification", "gridDesign accompanies explicit z; nPoints creates its own mode-root provenance.");
end

nAvailableModes = size(self.nativeModes,2);
if hasPointCount
    [z,maximumModeCount,gridDesign] = IMDiscreteTransformTools.pointsForExactCount(self,options.nPoints,nAvailableModes);
else
    z = options.z(:);
    maximumModeCount = min(nAvailableModes,length(z));
    gridDesign = IMDiscreteTransformTools.validatedGridDesign(self,z,options.gridDesign);
end

[searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = emptySearchColumns();
transform = [];
assessment = [];
for modeCount = maximumModeCount:-1:1
    try
        [candidateTransform,candidateAssessment] = self.discreteTransform(z=z,nModes=modeCount,variables=options.variables, ...
            gramTolerance=options.gramTolerance,allowRetainedPrefix=true);
        candidateGramError = candidateAssessment.gramPolicy.error(end);
        accepted = candidateAssessment.gramPolicy.accepted(end);
        [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
            appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
            "gramSearch",modeCount,true,accepted,candidateAssessment.retainedModeCount,candidateGramError,"","");
        if accepted
            transform = candidateTransform;
            assessment = candidateAssessment;
            break
        end
    catch exception
        if ~isRecoverableSearchFailure(exception)
            rethrow(exception)
        end
        [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
            appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
            "gramSearch",modeCount,false,false,0,NaN,string(exception.identifier),string(exception.message));
    end
end
if isempty(assessment)
    error("IMBasisSet:NoCertifiedDiscreteTransform", "No independently fitted family count passed the Gram tolerance on the supplied points.");
end

if ~isempty(options.leakageTolerance) || ~isempty(options.quadraticAliasingTolerance)
    modeCount = assessment.retainedModeCount;
    while true
        try
            [transform,assessment] = self.discreteTransform(z=z,nModes=modeCount,variables=options.variables, ...
                gramTolerance=options.gramTolerance,leakageTolerance=options.leakageTolerance, ...
                quadraticAliasingTolerance=options.quadraticAliasingTolerance,nCheckModes=options.nCheckModes, ...
                allowRetainedPrefix=true);
            accepted = assessment.retainedModeCount == modeCount;
            [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
                appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
                "policyRefit",modeCount,true,accepted,assessment.retainedModeCount,assessment.gramPolicy.error(end),"","");
            if accepted
                break
            end
            modeCount = assessment.retainedModeCount;
        catch exception
            if ~isRecoverableSearchFailure(exception) || modeCount == 1
                rethrow(exception)
            end
            [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
                appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
                "policyRefit",modeCount,false,false,0,NaN,string(exception.identifier),string(exception.message));
            modeCount = modeCount-1;
        end
    end
end

certificationSearch = table(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
    VariableNames=["stage","modeCount","fitSucceeded","accepted","retainedModeCount","gramError","failureIdentifier","failureMessage"]);
assessment = assessment.withCertificationMetadata(gridDesign,certificationSearch);
if assessment.weightFitModeCount ~= assessment.retainedModeCount
    error("IMBasisSet:InconsistentCertifiedTransform", "Certified weights must be fitted to the final retained family count.");
end
end

function [stage,modeCount,fitSucceeded,accepted,retainedModeCount,gramError,failureIdentifier,failureMessage] = emptySearchColumns()
stage = strings(0,1);
modeCount = zeros(0,1);
fitSucceeded = false(0,1);
accepted = false(0,1);
retainedModeCount = zeros(0,1);
gramError = zeros(0,1);
failureIdentifier = strings(0,1);
failureMessage = strings(0,1);
end

function [stage,modeCount,fitSucceeded,accepted,retainedModeCount,gramError,failureIdentifier,failureMessage] = appendSearch( ...
        stage,modeCount,fitSucceeded,accepted,retainedModeCount,gramError,failureIdentifier,failureMessage, ...
        newStage,newModeCount,newFitSucceeded,newAccepted,newRetainedModeCount,newGramError,newFailureIdentifier,newFailureMessage)
stage(end+1,1) = newStage;
modeCount(end+1,1) = newModeCount;
fitSucceeded(end+1,1) = newFitSucceeded;
accepted(end+1,1) = newAccepted;
retainedModeCount(end+1,1) = newRetainedModeCount;
gramError(end+1,1) = newGramError;
failureIdentifier(end+1,1) = newFailureIdentifier;
failureMessage(end+1,1) = newFailureMessage;
end

function tf = isRecoverableSearchFailure(exception)
tf = ismember(string(exception.identifier),["IMBasisSet:DiscreteTransformPolicyFailed", ...
    "IMBasisSet:NoAcceptableDiscreteTransformPrefix","IMBasisSet:QuadratureWeightFitFailed"]);
end
