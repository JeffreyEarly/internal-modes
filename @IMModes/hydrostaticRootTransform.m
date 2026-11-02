function [transform,grid,assessment] = hydrostaticRootTransform(background,nPoints,options)
% Build the rigid-lid hydrostatic F/G mode-root transform with a qualified band.
%
% The rigid-lid hydrostatic problem has non-active vertical boundaries,
%
% $$-G_j''=\frac{N^2}{g h_j}G_j,\qquad G_j(z_b)=G_j(z_s)=0,\qquad F_j=h_jG_j',$$
%
% normalized so $$\int (N^2/g)G_iG_j\,dz=\delta_{ij}$$ and therefore
% $$\int F_iF_j\,dz=h_j\delta_{ij}$$. The aligned family is
% $$F_0,\ldots,F_r$$ and $$G_0,\ldots,G_r$$, where the analytical barotropic
% sector $$F_0=1$$, $$G_0=0$$ is appended explicitly and $$G_0$$ is inactive.
%
% One G family is solved and reused. The $$N_z$$ sample points are both
% endpoints and the interior roots of $$G_{N_z-1}$$, as in Early, Lelong and
% Smith (2020). For each candidate retained count $$r\le N_z-2$$ one
% nonnegative full-depth weight vector $$w$$ is fitted to the stacked
% normalized F and G Gram errors of that band, with
% $$W_F=\mathrm{diag}(w)$$ and $$W_G=\mathrm{diag}(N^2(z_i)w_i/g)$$. The
% band is accepted when both
%
% $$\varepsilon_V=\left\|H_V^{-1/2}(\Phi_V^\mathsf{T}W_V\Phi_V-H_V)H_V^{-1/2}\right\|_2,\qquad V\in\{F,G\},$$
%
% are at most `gramTolerance` and both sampled Gram matrices have full rank.
% Here $$H_F=\mathrm{diag}(D,h_1,\ldots,h_r)$$ and the active $$H_G=I_r$$.
% $$\varepsilon_V$$ bounds the relative energy error of any combination of
% retained modes; the fitted weights are the rule for energies and other
% integrals on the grid.
%
% Coefficients are recovered by one of two forward operators, both exact for
% retained modes:
%
% - `"truncatedCollocation"` (default) inverts the square systems over the
%   whole root-grid span, $$F_0,\ldots,F_{N_z-1}$$ on all points and
%   $$G_1,\ldots,G_{N_z-2}$$ on the interior points, and keeps the rows
%   $$0,\ldots,r$$. Modes $$r+1,\ldots,N_z-1$$ (F) and $$r+1,\ldots,N_z-2$$ (G)
%   then produce no retained coefficients.
% - `"galerkin"` applies $$A_V=(\Phi_V^\mathsf{T}W_V\Phi_V)^{-1}\Phi_V^\mathsf{T}W_V$$,
%   computed by weighted QR, which lets those intermediate modes leak into
%   the retained band.
%
% Modes beyond the span alias under either operator. The search, its
% acceptance and the weights are the same for both.
%
% A specified `nInternalModes` is strict: that band is fitted and either
% accepted or rejected with diagnostics, never truncated. When it is
% omitted, bands are tested from $$N_z-2$$ downward and the largest accepted
% band is returned, with the search and the criterion that rejected the next
% larger band recorded in `assessment`. Acceptance is linear: it does not
% qualify nonlinear products, and it is not a claim of an optimal grid.
%
% Before fitting, the retained, omitted-check and generating modes are
% compared against an independent solve at `refinementFactor` times the
% resolution (`assessModeConvergence`); equivalent depths, F and G shapes
% and the root points must agree within `convergenceTolerance`. Without a
% supplied solver, the WKB spectral resolution starts at
% max(128,4*nPoints) and doubles until this check passes (at most 2048 or
% five attempts). A supplied solver is used as given.
%
% ```matlab
% background = IMModeBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
% [transform,grid,assessment] = IMModes.hydrostaticRootTransform(background,33);
% coefficients = transform.transformForward(F,variable="F");
% assessment.retainedInternalModeCount
% ```
%
% This method requires `lsqlin` from Optimization Toolbox for the weight fit.
%
% - Topic: Create modes
% - Declaration: [transform,grid,assessment] = IMModes.hydrostaticRootTransform(background,nPoints,options)
% - Parameter background: rigid-lid IMModeBackground with positive N2
% - Parameter nPoints: number of sample points, including both endpoints
% - Parameter options.nInternalModes: strict retained internal count r; empty selects automatically
% - Parameter options.solver: IMSolverSpectral used as given; default refines a WKB solver until converged
% - Parameter options.gramTolerance: acceptance tolerance for both F and G energy errors
% - Parameter options.convergenceTolerance: tolerance for the independent refinement check
% - Parameter options.refinementFactor: resolution ratio of the refinement solve
% - Parameter options.nOmittedModes: omitted internal modes used in sensitivity diagnostics
% - Parameter options.forwardOperator: `"truncatedCollocation"` (default) or `"galerkin"`
% - Returns transform: aligned IMInternalModesDiscreteTransform with columns F0..Fr
% - Returns grid: IMVerticalGrid holding the final points and fitted weights
% - Returns assessment: IMHydrostaticRootTransformAssessment
arguments (Input)
    background (1,1) IMModeBackground
    nPoints (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nPoints,3)}
    options.nInternalModes double {mustBeInteger,mustBePositive} = zeros(1,0)
    options.solver = []
    options.gramTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-2
    options.convergenceTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-6
    options.refinementFactor (1,1) double {mustBeReal,mustBeFinite,mustBeGreaterThan(options.refinementFactor,1)} = 1.5
    options.nOmittedModes (1,1) double {mustBeInteger,mustBeNonnegative} = 4
    options.forwardOperator (1,1) string {mustBeMember(options.forwardOperator,["truncatedCollocation","galerkin"])} = "truncatedCollocation"
end
arguments (Output)
    transform (1,1) IMInternalModesDiscreteTransform
    grid (1,1) IMVerticalGrid
    assessment (1,1) IMHydrostaticRootTransformAssessment
end

if background.surfaceCondition ~= "rigidLid"
    error("IMModes:HydrostaticRootRequiresRigidLid","The hydrostatic root transform has non-active boundaries; use a background with surfaceCondition=""rigidLid"".");
end
zDomain = background.zDomain;
D = background.D;
zCheck = linspace(zDomain(1),zDomain(2),2049).';
N2Check = background.N2(zCheck);
if any(~isfinite(N2Check)) || any(N2Check <= 0)
    error("IMModes:HydrostaticRootRequiresPositiveStratification","The hydrostatic root transform requires finite N2 > 0 throughout the domain.");
end
nMaximumInternal = nPoints-2;
strict = ~isempty(options.nInternalModes);
if strict
    if ~isscalar(options.nInternalModes) || options.nInternalModes > nMaximumInternal
        error("IMModes:InvalidHydrostaticRootBand","nInternalModes must be a scalar no larger than nPoints-2 = %d.",nMaximumInternal);
    end
    candidates = options.nInternalModes;
else
    candidates = nMaximumInternal:-1:1;
end

% Solve one G family: retained candidates, omitted checks and the generating
% mode. A default solver is refined automatically until an independent solve
% at refinementFactor times its resolution agrees; a supplied solver is used
% as given and must already pass that check.
args = namedargs2cell(struct("N2",background.N2,"zDomain",zDomain,"rotationRate",background.rotationRate,"latitude",background.latitude,"g",background.g));
evp = IMInternalModes.hydrostaticGModes(args{:});
nSolved = max(nPoints-1,nMaximumInternal+options.nOmittedModes);
normalization = "geostrophic";
if isempty(options.solver)
    resolutions = max(128,4*nPoints)*2.^(0:4);
    resolutions = resolutions(resolutions <= max(2048,max(128,4*nPoints)));
    coordinateKind = "wkb";
else
    if ~isa(options.solver,"IMSolverSpectral")
        error("IMModes:HydrostaticRootRequiresSpectralSolver","The hydrostatic root transform needs spectral root finding; supply an IMSolverSpectral.");
    end
    resolutions = options.solver.nEVP;
    coordinateKind = options.solver.coordinateKind;
end
attempts = zeros(0,4);
for nEVP = resolutions
    solver = IMSolverSpectral(nEVP=nEVP,coordinateKind=coordinateKind);
    refinedSolver = IMSolverSpectral(nEVP=ceil(options.refinementFactor*nEVP),coordinateKind=coordinateKind);
    [basis,convergence,z,rootDesign] = convergedFamily(evp,solver,refinedSolver,nSolved,nPoints,normalization,zDomain);
    attempts(end+1,:) = [nEVP convergence.maximumEquivalentDepthError convergence.maximumShapeError convergence.rootDisplacement]; %#ok<AGROW>
    converged = convergence.maximumEquivalentDepthError <= options.convergenceTolerance ...
        && convergence.maximumShapeError <= options.convergenceTolerance && convergence.rootDisplacement <= options.convergenceTolerance;
    if converged
        break
    end
end
convergence.tolerance = options.convergenceTolerance;
convergence.attempts = array2table(attempts,VariableNames=["nEVP","maximumEquivalentDepthError","maximumShapeError","rootDisplacement"]);
if ~converged
    error("IMModes:UnconvergedHydrostaticRootModes", ...
        "Modes 1-%d are not converged between nEVP=%d and %d (equivalent depth %.3g, shape %.3g, root displacement %.3g; tolerance %.3g). Supply an IMSolverSpectral with larger nEVP.", ...
        nSolved,convergence.candidateNEVP,convergence.referenceNEVP,convergence.maximumEquivalentDepthError,convergence.maximumShapeError,convergence.rootDisplacement,options.convergenceTolerance);
end
z([1 end]) = zDomain;

% Sampled aligned family with the analytical barotropic sector.
h = basis.h(1:nSolved);
FSampled = [ones(nPoints,1) basis.F(z,normalization=normalization)];
GSampled = [zeros(nPoints,1) basis.G(z,normalization=normalization)];
GSampled([1 end],:) = 0;
HF = [D h(:).'];
HG = [0 ones(1,nSolved)];
N2Values = background.N2(z);
GMetric = N2Values(:)/background.g;
geometricWeights = diff([zDomain(1);0.5*(z(1:end-1)+z(2:end));zDomain(2)]);
% Square collocation over the full root-grid span: F0..F_{Nz-1} on all
% points and G_1..G_{Nz-2} on the interior points, where G vanishes at both
% endpoints. Rows are truncated to the retained band afterwards.
[collocationF,collocationG,collocationDiagnostics] = collocationInverses(FSampled(:,1:nPoints),GSampled(2:end-1,2:nPoints-1));

searchRows = cell(numel(candidates),1);
accepted = [];
fits = cell(numel(candidates),1);
for iCandidate = 1:numel(candidates)
    r = candidates(iCandidate);
    [w,fitResidual] = fitBandWeights(FSampled(:,1:r+1),GSampled(:,2:r+1),HF(1:r+1),GMetric,geometricWeights,D);
    row = bandDiagnostics(r,w,fitResidual,FSampled,GSampled,HF,GMetric,options.gramTolerance,options.nOmittedModes,D, ...
        options.forwardOperator,collocationF,collocationG);
    searchRows{iCandidate} = row;
    fits{iCandidate} = w;
    if row.accepted
        accepted = iCandidate;
        break
    end
end
search = struct2table([searchRows{~cellfun(@isempty,searchRows)}]);
if isempty(accepted)
    last = search(end,:);
    if strict
        error("IMModes:HydrostaticRootBandRejected", ...
            "The requested band r=%d fails: F energy error %.3g, G energy error %.3g (tolerance %.3g), F rank %d/%d, G rank %d/%d. Request fewer modes or omit nInternalModes for automatic selection.", ...
            last.internalModeCount,last.energyErrorF,last.energyErrorG,options.gramTolerance,last.rankF,last.internalModeCount+1,last.rankG,last.internalModeCount);
    end
    error("IMModes:NoHydrostaticRootBand","No retained band passed the F and G energy tolerance %.3g on the %d-point root grid.",options.gramTolerance,nPoints);
end
r = candidates(accepted);
w = fits{accepted};
if strict
    selection = "strict";
    limitingCriterion = "requested";
elseif accepted == 1
    selection = "automatic";
    limitingCriterion = "maximumRootGridBand";
else
    selection = "automatic";
    limitingCriterion = search.rejectionReason(accepted-1);
end

columns = 1:r+1;
modeNumber = 0:r;
FColumns = FSampled(:,columns);
GColumns = GSampled(:,columns);
channel = @(metric,target,active) struct(available=true,reason="",activeModeMask=active, ...
    metricMatrix=diag(metric),targetGramMatrix=diag(target),targetMajorantGramMatrix=diag(target));
channelData = struct(F=channel(w,HF(columns),true(1,r+1)),G=channel(GMetric.*w,HG(columns),[false true(1,r)]));
if options.forwardOperator == "truncatedCollocation"
    [forwardF,forwardG] = truncatedCollocation(r,collocationF,collocationG);
    channelData.F.suppliedForwardMatrix = forwardF;
    channelData.F.suppliedForwardOperator = "truncatedCollocation";
    channelData.G.suppliedForwardMatrix = forwardG;
    channelData.G.suppliedForwardOperator = "truncatedCollocation";
end
metadata = struct(construction="IMModes.hydrostaticRootTransform",background=background,surfaceCondition="rigidLid", ...
    columnLabels=["barotropic:0" "internal:"+string(1:r)],relation="F = h dG/dz",generatingMode=nPoints-1);
transform = IMInternalModesDiscreteTransform(z=z,weights=w,modeNumber=modeNumber,h=[Inf h(1:r)],normalization=normalization, ...
    inverseF=FColumns,inverseG=GColumns,endpointF=FColumns([end 1],:),endpointG=GColumns([end 1],:),channelData=channelData, ...
    primaryVariable="G",zDomain=zDomain,g=background.g,modeFamily="hydrostatic",N2Values=N2Values,problemMetadata=metadata);
gridDesign = struct(rule="hydrostaticRootTransform",generatingMode=nPoints-1,rootDesign=rootDesign, ...
    retainedInternalModeCount=r,weightRule="Nonnegative full-depth fit of stacked normalized F/G Gram errors for the retained band");
grid = IMVerticalGrid(z=z,weights=w,zDomain=zDomain,kind="modeRoots",design=gridDesign);
assessment = IMHydrostaticRootTransformAssessment(selection=selection,nPoints=nPoints,retainedInternalModeCount=r, ...
    gramTolerance=options.gramTolerance,limitingCriterion=limitingCriterion,search=search,convergence=convergence, ...
    forwardOperator=options.forwardOperator,collocation=collocationDiagnostics, ...
    final=search(accepted,:));
end

function [basis,convergence,z,rootDesign] = convergedFamily(evp,solver,refinedSolver,nSolved,nPoints,normalization,zDomain)
basis = solver.solveEVP(evp,nModes=nSolved);
refined = refinedSolver.solveEVP(evp,nModes=nSolved);
[zFine,wFine] = physicalReferenceQuadrature(zDomain);
labels = "internal:"+string(1:nSolved);
identity = struct(family="hydrostaticRootG",columnLabels=labels,normalization=normalization,zDomain=zDomain);
candidateModes = struct(identity=identity,values=struct(F=basis.F(zFine,normalization=normalization),G=basis.G(zFine,normalization=normalization)), ...
    equivalentDepths=basis.h(1:nSolved),provenance=struct(nEVP=solver.nEVP,coordinateKind=solver.coordinateKind));
referenceModes = struct(identity=identity,values=struct(F=refined.F(zFine,normalization=normalization),G=refined.G(zFine,normalization=normalization)), ...
    equivalentDepths=refined.h(1:nSolved),provenance=struct(nEVP=refinedSolver.nEVP,coordinateKind=refinedSolver.coordinateKind));
measurements = assessModeConvergence(candidateModes,referenceModes,zFine,wFine).measurements;
[z,rootDesign] = basis.pointsFromModeRoots(nModes=nPoints-2);
zRefined = refined.pointsFromModeRoots(nModes=nPoints-2);
if numel(z) ~= nPoints || numel(zRefined) ~= nPoints
    error("IMModes:HydrostaticRootCount","The generating mode G_%d must give exactly %d points including both endpoints; found %d.",nPoints-1,nPoints,numel(z));
end
convergence = struct(candidateNEVP=solver.nEVP,referenceNEVP=refinedSolver.nEVP,coordinateKind=solver.coordinateKind, ...
    maximumEquivalentDepthError=max(measurements.value(measurements.quantity=="equivalentDepth")), ...
    maximumShapeError=max(measurements.value(measurements.quantity=="shape")), ...
    rootDisplacement=max(abs(z-zRefined))/diff(zDomain),measurements=measurements);
end

function [z,w] = physicalReferenceQuadrature(zDomain)
grid = IMVerticalGrid.physicalDepth(zDomain,1025);
z = grid.z; w = grid.weights;
end

function [w,residual] = fitBandWeights(F,G,HF,GMetric,geometricWeights,depth)
% Stacked normalized Gram rows; off-diagonal rows carry sqrt(2) so the
% least-squares objective equals the full normalized Frobenius error.
[AF,bF] = gramRows(F,HF,ones(size(GMetric)));
[AG,bG] = gramRows(G,ones(1,size(G,2)),GMetric);
A = [AF;AG]; b = [bF;bG];
w = IMDiscreteTransformTools.fitQuadrature(A,b,geometricWeights,depth,true,true);
residual = norm(A*w-b);
end

function [A,b] = gramRows(Phi,H,metric)
m = size(Phi,2);
[I,J] = find(triu(true(m)));
scale = 1./sqrt(H(I).*H(J));
offDiagonal = I ~= J;
scale(offDiagonal) = sqrt(2)*scale(offDiagonal);
A = (Phi(:,I).*Phi(:,J).*metric).'.*scale(:);
b = double(I == J);
end

function row = bandDiagnostics(r,w,fitResidual,F,G,HF,GMetric,tolerance,nOmitted,depth,forwardOperator,collocationF,collocationG)
retainedF = 1:r+1;
retainedG = 2:r+1;
[epsF,rankF,condF,AF] = channelError(F(:,retainedF),w,HF(retainedF));
[epsG,rankG,condG,AG] = channelError(G(:,retainedG),GMetric.*w,ones(1,r));
if forwardOperator == "truncatedCollocation"
    [AF,AGFull] = truncatedCollocation(r,collocationF,collocationG);
    AG = AGFull(2:end,:);
end
omitted = r+2:min(size(F,2),r+1+nOmitted);
if isempty(omitted)
    gainF = NaN; gainG = NaN;
else
    gainF = norm(sqrt(HF(retainedF)).'.*(AF*F(:,omitted))./sqrt(HF(omitted)),2);
    gainG = norm(AG*G(:,omitted),2);
end
modalMean = max([0 abs(w.'*F(:,2:r+1))./sqrt(depth*HF(2:r+1))]);
fullRank = rankF == r+1 && rankG == r;
acceptedBand = fullRank && epsF <= tolerance && epsG <= tolerance;
reason = "";
if ~acceptedBand
    if rankF < r+1, reason = "rankF";
    elseif rankG < r, reason = "rankG";
    elseif epsF > tolerance && epsG > tolerance, reason = "energyErrorFG";
    elseif epsF > tolerance, reason = "energyErrorF";
    else, reason = "energyErrorG";
    end
end
row = struct(internalModeCount=r,energyErrorF=epsF,energyErrorG=epsG,rankF=rankF,rankG=rankG, ...
    gramConditionF=condF,gramConditionG=condG,minimumWeight=min(w),maximumWeight=max(w), ...
    zeroWeightCount=nnz(w <= 1e-12*depth),maximumModalMean=modalMean,omittedModeCount=numel(omitted), ...
    omittedGainF=gainF,omittedGainG=gainG,fitResidual=fitResidual,accepted=acceptedBand,rejectionReason=reason);
end

function [inverseF,inverseG,diagnostics] = collocationInverses(F,G)
% Column-equilibrated square inverses, without forming an explicit inverse
% of the unscaled matrices.
scaleF = max(abs(F),[],1);
scaleG = max(abs(G),[],1);
rcondF = rcond(F./scaleF);
rcondG = rcond(G./scaleG);
if ~(rcondF > 1e3*eps) || ~(rcondG > 1e3*eps)
    error("IMModes:SingularHydrostaticRootCollocation","The square root-grid collocation system is numerically singular (reciprocal condition F %.3g, G %.3g).",rcondF,rcondG);
end
inverseF = ((F./scaleF)\eye(size(F,1)))./scaleF.';
inverseG = ((G./scaleG)\eye(size(G,1)))./scaleG.';
diagnostics = struct(conditionF=cond(F./scaleF),conditionG=cond(G./scaleG),columnScaling="maximum absolute value", ...
    spanF="F_0..F_{Nz-1} on all points",spanG="G_1..G_{Nz-2} on interior points");
end

function [forwardF,forwardG] = truncatedCollocation(r,inverseF,inverseG)
nPoints = size(inverseF,2);
forwardF = inverseF(1:r+1,:);
forwardG = zeros(r+1,nPoints);
forwardG(2:r+1,2:end-1) = inverseG(1:r,:);
end

function [epsilon,rankValue,conditionNumber,forward] = channelError(Phi,metric,H)
Gamma = Phi.'*(metric.*Phi);
Gamma = 0.5*(Gamma+Gamma.');
scale = 1./sqrt(H(:));
epsilon = norm(scale.*(Gamma-diag(H)).*scale.',2);
singularValues = svd(Gamma);
rankValue = sum(singularValues > numel(H)*eps(max(singularValues)));
conditionNumber = singularValues(1)/singularValues(end);
rootMetric = sqrt(metric);
[Q,R] = qr(rootMetric.*Phi,0);
forward = R\(Q.'.*rootMetric.');
end
