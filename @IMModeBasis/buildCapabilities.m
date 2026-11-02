function description = buildCapabilities(self,scientificInput)
% Describe supported continuous operations without preparing any operator.
%
% The immutable description is shared by scalar and collection validation.
% Normalization names describe structurally available rules; each retained
% column must still have a finite nonzero magnitude in the selected rule.
% Recovery availability is a property of a prepared observation system,
% rather than of a continuous profile or physical pairing.
%
% - Topic: Developer topics
% - Returns description: continuous capabilities and conditional preparation requirements
% - Developer: true
variables = ["F","G","APV"];
variables = variables(ismember(variables,self.profileVariables_));
orders = repmat({0},1,numel(variables));
for iVariable = 1:numel(variables)
    if variables(iVariable) == "F" && (~isempty(self.source) || self.waveConvention || ismember("Fz",self.profileVariables_))
        orders{iVariable} = 0:1;
    elseif variables(iVariable) == "G"
        orders{iVariable} = find(ismember(["G","Gz","Gzz"],self.profileVariables_))-1;
    end
end
aliases = struct("variable",{},"baseVariable",{},"derivativeOffset",{});
if ismember("F",variables)
    aliases = struct(variable="Fz",baseVariable="F",derivativeOffset=1);
end
if ismember("G",variables)
    aliases = [aliases,struct("variable",{"Gz","Gzz"},"baseVariable",{"G","G"},"derivativeOffset",{1,2})];
end
isIntervalBasis = isa(self.source,"IMNeutralLabelModesBasis") || ~isempty(self.integrationIntervals_);
if isIntervalBasis
    derivativeConvention = "classicalIntervalInterior";
    traceConvention = "insideIntervalEndpoints";
else
    derivativeConvention = "classicalPhysical";
    traceConvention = "physicalDomainEvaluation";
end
spectrumNames = strings(1,0);
if self.hasEquivalentDepths_, spectrumNames(end+1) = "h"; end
if self.hasEigenvalues_, spectrumNames(end+1) = "eigenvalues"; end
if self.horizontalKind == "perColumn", spectrumNames(end+1) = "kSquared"; end
if self.hasFrequency_, spectrumNames(end+1) = "frequency"; end
evp = [];
if isfield(scientificInput,"evp")
    evp = scientificInput.evp;
elseif ~isempty(self.source) && isprop(self.source,"evp")
    evp = self.source.evp;
end
normalizationNames = strings(1,0);
if ismember("F",variables), normalizationNames = [normalizationNames,"depthRMSF","maxAbsF"]; end
if ismember("G",variables), normalizationNames = [normalizationNames,"depthRMSG","maxAbsG"]; end
if self.waveConvention && self.horizontalKind == "fixedWavenumber", normalizationNames(end+1) = "waveEnergy"; end
hasGeneralizedEnergy = ~isempty(self.sourceUnityFactors_) && ismethod(self.source,"gramMatrix");
if ~isempty(evp), hasGeneralizedEnergy = evp.innerProduct("G").hasInnerProduct; end
if hasGeneralizedEnergy && ~self.waveConvention
    normalizationNames(end+1) = "generalizedEnergy";
end
hasGeneralizedEnstrophy = ~isempty(self.sourceUnityFactors_) && ismember("generalizedPotentialEnstrophy",self.source.normalizationNames());
if ~isempty(evp), hasGeneralizedEnstrophy = ismember("generalizedPotentialEnstrophy",IMInternalModesNormalizationTools.available(evp)); end
if hasGeneralizedEnstrophy
    normalizationNames(end+1) = "generalizedPotentialEnstrophy";
end
if isa(self.source,"IMGeostrophicZeroAPVModesBasis") || isa(self.source,"IMAnalyticalGeostrophicZeroAPVModesBasis") || self.family == "geostrophicZeroAPV"
    normalizationNames(end+1) = "endpointResponse";
end
normalizationNames = unique(normalizationNames,"stable");
operations = struct("variable",{},"kind",{},"availability",{},"requirements",{});
if self.family == "hydrostaticF" && self.background.surfaceCondition == "rigidLid"
    for variable = ["F","G"]
        operations(end+1) = struct(variable=variable,kind="sampledGalerkin",availability="conditional",requirements=["fixed physical samples and dz quadrature","continuous physical pairing","separate prepared-recovery rank assessment"]); %#ok<AGROW>
    end
elseif isa(self.source,"IMNeutralLabelModesBasis")
    operations = struct(variable="APV",kind="prescribedContinuousDual",availability="conditional",requirements=["APV observations with fixed dz quadrature","nonzero f0 and representable continuous APV norms for retained columns","empty inventories preserve an empty physical dual","particle-label information for a complete neutral-label state"]);
end
recovery = struct(synthesisVariables=variables,preparation="IMDiscreteTransform.fromObservationMatrices",availability="requiresPreparedSystem",requirements=["explicit observation model and column identities","explicit scientific column activity","declared fitting or observation factor, or physical dual","successful rank and overlap assessment of the prepared system"],noiseEvidence="requiresDeclaredObservationFactor");
description = struct(variables=variables,derivativeOrders={orders},aliases=aliases,derivativeConvention=derivativeConvention,traceConvention=traceConvention,traceOrders={orders},spectrumNames=spectrumNames,normalizationNames=normalizationNames,normalizationRequirements="Every retained column requires a finite positive magnitude; signed rules preserve a separate norm signature.",projectionOperations=operations,sampledRecovery=recovery,supportsEmptyColumns=true);
end
