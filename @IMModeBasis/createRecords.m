function [scientific,construction] = createRecords(self,scientificInput,constructionInput)
% Snapshot the continuous specification separately from its construction.
%
% These value records describe the actual source and retained columns.
% Construction diagnostics do not qualify convergence, observation grids,
% or a numerical recovery operator. Requested settings remain available
% after selecting or reordering the returned columns.
arguments (Input)
    self (1,1) IMModeBasis
    scientificInput (1,1) struct
    constructionInput (1,1) struct
end
arguments (Output)
    scientific (1,1) struct
    construction (1,1) struct
end
background = self.background;
physical = struct(N2=background.N2,zDomain=background.zDomain,rotationRate=background.rotationRate,latitude=background.latitude,f0=background.f0,g=background.g,surfaceCondition=background.surfaceCondition,profileKind=background.profileKind,profileParameters=background.profileParameters);
requested = struct();
if isfield(scientificInput,"requested"), requested = scientificInput.requested; end
evp = [];
if isfield(scientificInput,"evp")
    evp = scientificInput.evp;
elseif ~isempty(self.source) && isprop(self.source,"evp")
    evp = self.source.evp;
end
evpRecord = struct();
pairings = struct();
if ~isempty(evp)
    evpRecord = struct(className=string(class(evp)),name=evp.name,formulation=evp.formulation,modeFamily=evp.modeFamily,parameters=evp.parameters,p=evp.p,q=evp.q,r=evp.r,surfaceBoundary=boundaryRecord(evp.surfaceBoundary),bottomBoundary=boundaryRecord(evp.bottomBoundary));
    pairings.F = evp.innerProduct("F");
    pairings.G = evp.innerProduct("G");
end
energyBoundary = struct();
energyBoundaryRole = "notSpecified";
if isfield(scientificInput,"energyBoundary")
    boundary = scientificInput.energyBoundary;
    energyBoundary = struct(g0=boundary.g0,gd=boundary.gd,activeEndpoints=boundary.activeEndpoints);
    energyBoundaryRole = "definesPairingAndEndpointCoordinates";
elseif ~isempty(evp) && isfield(evp.parameters,"g0")
    energyBoundary = struct(g0=evp.parameters.g0,gd=evp.parameters.gd);
    energyBoundaryRole = "definesPairingAndEndpointCoordinates";
elseif isa(self.source,"IMNeutralLabelModesBasis")
    energyBoundary = struct(g0=self.source.g0,gd=self.source.gd);
    energyBoundaryRole = "provenanceOnly";
end
intervals = self.integrationIntervals_;
intervalKind = "none";
if isfield(scientificInput,"intervals"), intervals = scientificInput.intervals; end
if isfield(scientificInput,"intervalKind"), intervalKind = scientificInput.intervalKind; end
intervalBoundaryTypes = strings(0,2);
endpointDerivativeTrace = "notApplicable";
if isa(self.source,"IMNeutralLabelModesBasis")
    intervals = self.source.neutralIntervals;
    intervalKind = "neutral";
    energyBoundaryRole = "provenanceOnly";
end
if ismember(intervalKind,["neutral","critical"])
    intervalBoundaryTypes = repmat("Dirichlet",size(intervals));
    endpointDerivativeTrace = "insideDeclaredInterval";
    pairings.G = struct(kind="depthL2",interiorWeight=1/background.D,hasInnerProduct=true,physicalEnergy=false);
    if intervalKind == "neutral"
        pairings.APV = struct(kind="prescribedContinuousAPVDual",profileRelation="APV=-f0*Gz",coefficientNorm="integral APV_j^2 dz",hasInnerProduct=true,isDegenerate=background.f0==0,physicalEnergy=false);
    end
end
endpoints = strings(1,0);
if isfield(scientificInput,"endpoints")
    endpoints = scientificInput.endpoints;
elseif ~isempty(self.source) && isprop(self.source,"endpoints")
    endpoints = self.source.endpoints;
end
if self.family == "geostrophicZeroAPV"
    energyBoundaryRole = "selectsEndpointCoordinatesOnly";
    pairings.endpointResponse = struct(kind="canonicalEndpointCoordinates",endpoints=endpoints,normalization="unitEndpointResponse");
end
profileRelation = "sourceDefined";
if self.waveConvention
    profileRelation = "F=D*Gz";
elseif self.family == "meanDensityAnomaly"
    profileRelation = "F(z)=integral_z^surface N2*G/g dz";
elseif self.family == "neutralLabel"
    profileRelation = "APV=-f0*Gz";
elseif self.family == "geostrophicZeroAPV"
    profileRelation = "Gz=-g*k^2*F/f0^2; Fz=-N2*G/g";
elseif ismember(self.family,["hydrostaticF","hydrostaticG","geostrophicAPV","geostrophicGeneralizedPotentialEnstrophy"])
    profileRelation = "Fz=-N2*G/g";
end
profileBoundarySemantics = "sourceDefined";
if ~isempty(evp), profileBoundarySemantics = "canonicalSourceBoundaryConditions"; end
if ismember(intervalKind,["neutral","critical"])
    profileBoundarySemantics = "G=0 at both endpoints of each declared interval; zero extension outside";
elseif self.family == "geostrophicZeroAPV"
    profileBoundarySemantics = "canonical surface/bottom endpoint responses in the physical surface convention";
elseif self.family == "hydrostaticF" && isempty(self.source) && self.nColumns == 1 && self.modeInfo.branch == "barotropic"
    profileBoundarySemantics = "analytical rigid-lid F0=1, G0=0 before normalization";
    pairings.F = struct(kind="L2",interiorWeight=1,hasInnerProduct=true);
    pairings.G = struct(kind="hydrostaticG",interiorWeight=@(z) background.N2(z)/background.g,hasInnerProduct=true);
end
normalization = struct(name=self.normalization,scales=self.scaleFactors_,normSignatures=self.normSignatures);
horizontal = struct(kind=self.horizontalKind,k=self.k,kSquared=self.kSquared);
includesAnalyticalBarotropic = self.family == "hydrostaticF" && any(self.modeInfo.branch == "barotropic");
resolved = struct(evp=evpRecord,pairings=pairings,energyBoundary=energyBoundary,energyBoundaryRole=energyBoundaryRole,intervals=intervals,intervalKind=intervalKind,intervalBoundaryTypes=intervalBoundaryTypes,endpointDerivativeTrace=endpointDerivativeTrace,endpoints=endpoints,normalization=normalization,profileRelation=profileRelation,profileBoundarySemantics=profileBoundarySemantics,horizontal=horizontal,columnLabels=self.columnLabels,modeInfo=self.modeInfo,sourceColumns=self.columns,includesAnalyticalBarotropic=includesAnalyticalBarotropic);
scientific = struct(schemaVersion=1,family=self.family,background=physical,requested=requested,resolved=resolved);

requestedSolver = struct();
if isfield(constructionInput,"requestedSolver")
    requestedSolver = solverRecord(constructionInput.requestedSolver,false);
end
actualSolver = struct();
if ~isempty(self.source) && isprop(self.source,"solver")
    actualSolver = solverRecord(self.source.solver,true);
elseif isfield(constructionInput,"actualSolver")
    actualSolver = solverRecord(constructionInput.actualSolver,false);
end
backendKind = "intrinsic";
if isfield(constructionInput,"backendKind")
    backendKind = constructionInput.backendKind;
elseif isa(self.source,"IMBasisSet") || isa(self.source,"IMGeostrophicZeroAPVModesBasis")
    backendKind = "numerical";
elseif ~isempty(self.source)
    backendKind = "analytical";
end
sourceClass = "";
if ~isempty(self.source), sourceClass = string(class(self.source)); end
if isfield(constructionInput,"sourceClass"), sourceClass = constructionInput.sourceClass; end
constructionSkipped = false;
if isfield(constructionInput,"constructionSkipped"), constructionSkipped = constructionInput.constructionSkipped; end
resolvedConstruction = struct(backendKind=backendKind,sourceClass=sourceClass,solver=actualSolver,constructionSkipped=constructionSkipped);
selection = struct(request=struct(),sourceEigenpairIndices=zeros(1,0),sourceModeNumbers=zeros(1,0),diagnostics=struct(),activeSourceColumns=self.columns,columnLabels=self.columnLabels,initialColumnIndices=1:self.nColumns);
if isfield(constructionInput,"selection")
    raw = constructionInput.selection;
    selection.sourceEigenpairIndices = reshape(raw.sortIndex,1,[]);
    selection.sourceModeNumbers = reshape(raw.modeNumber,1,[]);
    selection.diagnostics = raw.modeSelectionDiagnostics;
end
if isfield(constructionInput,"request"), selection.request = constructionInput.request; end
sourceDiagnostics = struct();
if ~isempty(self.source) && isprop(self.source,"metadata"), sourceDiagnostics = self.source.metadata; end
evidence = struct(solveAccuracy=struct(status="unverified",reason="No independent solve-accuracy qualification is attached."),referenceConvergence=struct(status="unverified",reason="Construction and retained counts do not establish convergence."),sourceDiagnostics=sourceDiagnostics);
construction = struct(schemaVersion=1,requested=struct(solver=requestedSolver),resolved=resolvedConstruction,selection=selection,evidence=evidence);
end

function record = boundaryRecord(boundary)
record = struct(a=boundary.a,b=boundary.b,c=boundary.c,d=boundary.d);
end

function record = solverRecord(solver,hasNativeRepresentation)
record = struct(className=string(class(solver)),nEVP=[],coordinateKind="notApplicable",zNative=[],xNative=[],nativeRepresentationAvailable=false);
if isprop(solver,"nEVP"), record.nEVP = solver.nEVP; end
if isprop(solver,"coordinateKind")
    record.coordinateKind = solver.coordinateKind;
elseif isa(solver,"IMSolverFiniteDifference")
    record.coordinateKind = "z";
elseif isa(solver,"IMSolverAnalytical")
    record.coordinateKind = "analytical";
end
% A finite-difference grid is the requested representation itself.
hasNativeRepresentation = hasNativeRepresentation || isa(solver,"IMSolverFiniteDifference");
if hasNativeRepresentation
    if isprop(solver,"zNative"), record.zNative = solver.zNative; end
    if isprop(solver,"xNative"), record.xNative = solver.xNative; end
    record.nativeRepresentationAvailable = ~isempty(record.zNative);
end
end
