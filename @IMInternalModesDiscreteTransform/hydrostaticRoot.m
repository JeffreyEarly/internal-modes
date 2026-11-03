function [transform,grid,assessment,basis] = hydrostaticRoot(background,nPoints,options)
% Build the rigid-lid hydrostatic F/G mode-root transform with a qualified band.
%
% The rigid-lid hydrostatic problem has non-active vertical boundaries,
%
% $$-G_j''=\frac{N^2}{g h_j}G_j,\qquad G_j(z_b)=G_j(z_s)=0,\qquad F_j=h_jG_j',$$
%
% using the shared `depthRMSF` normalization:
% $$\int F_iF_j\,dz=D\delta_{ij}$$ and
% $$\int (N^2/g)G_iG_j\,dz=(D/h_j)\delta_{ij}$$. The aligned family is
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
% Here $$H_F=D I_{r+1}$$ and the active $$H_G=\mathrm{diag}(D/h_1,\ldots,D/h_r)$$.
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
% five resolution levels), trying physical z as a fallback at each level.
% Root-count failures are recorded and refined. A supplied solver is used as given.
%
% ```matlab
% background = IMBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
% [transform,grid,assessment] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33);
% coefficients = transform.transformForward(F,variable="F");
% assessment.retainedInternalModeCount
% ```
%
% This method requires `lsqlin` from Optimization Toolbox for the weight fit.
%
% - Topic: Create discrete transforms
% - Declaration: [transform,grid,assessment,basis] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,nPoints,options)
% - Parameter background: rigid-lid IMBackground with positive N2
% - Parameter nPoints: number of sample points, including both endpoints
% - Parameter options.nInternalModes: strict retained internal count r; empty selects automatically
% - Parameter options.solver: IMSolverSpectral used as given; default refines WKB with physical-z fallback
% - Parameter options.gramTolerance: acceptance tolerance for both F and G energy errors
% - Parameter options.convergenceTolerance: tolerance for the independent refinement check
% - Parameter options.refinementFactor: resolution ratio of the refinement solve
% - Parameter options.nOmittedModes: omitted internal modes used in sensitivity diagnostics
% - Parameter options.forwardOperator: `"truncatedCollocation"` (default) or `"galerkin"`
% - Returns transform: aligned IMInternalModesDiscreteTransform with columns F0..Fr
% - Returns grid: IMVerticalGrid holding the final points and fitted weights
% - Returns assessment: IMHydrostaticRootTransformAssessment
% - Returns basis: retained native IMHydrostaticModesBasis, also available as transform.basis
arguments (Input)
    background (1,1) IMBackground
    nPoints (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(nPoints,3)}
    options.nInternalModes double {mustBeInteger,mustBeNonnegative,mustBeFinite} = zeros(1,0)
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
    basis (1,1) IMHydrostaticModesBasis
end

if background.surfaceCondition ~= "rigidLid"
    error("IMInternalModesDiscreteTransform:HydrostaticRootRequiresRigidLid","The hydrostatic root transform has non-active boundaries; use a background with surfaceCondition=""rigidLid"".");
end
zDomain = background.zDomain;
D = background.D;
zCheck = linspace(zDomain(1),zDomain(2),2049).';
N2Check = background.N2(zCheck);
if any(~isfinite(N2Check)) || any(N2Check <= 0)
    error("IMInternalModesDiscreteTransform:HydrostaticRootRequiresPositiveStratification","The hydrostatic root transform requires finite N2 > 0 throughout the domain.");
end
nMaximumInternal = nPoints-2;
strict = ~isempty(options.nInternalModes);
if strict
    if ~isscalar(options.nInternalModes) || options.nInternalModes > nMaximumInternal
        error("IMInternalModesDiscreteTransform:InvalidHydrostaticRootBand","nInternalModes must be a scalar no larger than nPoints-2 = %d.",nMaximumInternal);
    end
    candidates = options.nInternalModes;
else
    candidates = nMaximumInternal:-1:0;
end

% Solve one G family: retained candidates, omitted checks and the generating
% mode. A default solver is refined automatically until an independent solve
% at refinementFactor times its resolution agrees; a supplied solver is used
% as given and must already pass that check.
nSolved = max(nPoints-1,nMaximumInternal+options.nOmittedModes);
if isempty(options.solver)
    resolutions = max(128,4*nPoints)*2.^(0:4);
    resolutions = resolutions(resolutions <= max(2048,max(128,4*nPoints)));
    coordinateKinds = ["wkb","z"];
else
    if ~isa(options.solver,"IMSolverSpectral")
        error("IMInternalModesDiscreteTransform:HydrostaticRootRequiresSpectralSolver","The hydrostatic root transform needs spectral root finding; supply an IMSolverSpectral.");
    end
    resolutions = options.solver.nEVP;
    coordinateKinds = string(options.solver.coordinateKind);
end
attempts = struct([]);
converged = false;
for nEVP = resolutions
    for coordinateKind = coordinateKinds
        if isempty(options.solver)
            solver = IMSolverSpectral(nEVP=nEVP,coordinateKind=coordinateKind);
        else
            solver = options.solver;
        end
        refinedSolver = IMSolverSpectral(nEVP=ceil(options.refinementFactor*nEVP),coordinateKind=coordinateKind);
        try
            [solvedBasis,convergence,z,rootDesign] = convergedFamily(background,solver,refinedSolver,nSolved,nPoints,zDomain);
            attempt = struct(nEVP=nEVP,coordinateKind=coordinateKind,maximumEquivalentDepthError=convergence.maximumEquivalentDepthError, ...
                maximumShapeError=convergence.maximumShapeError,rootDisplacement=convergence.rootDisplacement,status="measured",errorIdentifier="");
            converged = convergence.maximumEquivalentDepthError <= options.convergenceTolerance ...
                && convergence.maximumShapeError <= options.convergenceTolerance && convergence.rootDisplacement <= options.convergenceTolerance;
        catch exception
            recoverable = ismember(string(exception.identifier),["IMInternalModesDiscreteTransform:HydrostaticRootCount", ...
                "IMBasisSet:InsufficientQuadraturePoints","IMSolver:InsufficientModes"]);
            if ~recoverable, rethrow(exception); end
            attempt = struct(nEVP=nEVP,coordinateKind=coordinateKind,maximumEquivalentDepthError=Inf, ...
                maximumShapeError=Inf,rootDisplacement=Inf,status="underresolved",errorIdentifier=string(exception.identifier));
            convergence = struct(candidateNEVP=nEVP,referenceNEVP=refinedSolver.nEVP,coordinateKind=coordinateKind, ...
                maximumEquivalentDepthError=Inf,maximumShapeError=Inf,rootDisplacement=Inf,measurements=table());
        end
        attempts = [attempts;attempt]; %#ok<AGROW>
        if converged, break, end
    end
    if converged, break, end
end
convergence.tolerance = options.convergenceTolerance;
convergence.attempts = struct2table(attempts);
if ~converged
    error("IMInternalModesDiscreteTransform:UnconvergedHydrostaticRootModes", ...
        "Modes 1-%d or generating G_%d roots could not be qualified after %d attempts; last coordinate %s, nEVP=%d/%d (depth %.3g, shape %.3g, root displacement %.3g; tolerance %.3g). Supply IMSolverSpectral with a suitable coordinateKind and larger nEVP. See coordinateKind=""z"" for strongly stretched backgrounds.", ...
        nSolved,nPoints-1,numel(attempts),convergence.coordinateKind,convergence.candidateNEVP,convergence.referenceNEVP, ...
        convergence.maximumEquivalentDepthError,convergence.maximumShapeError,convergence.rootDisplacement,options.convergenceTolerance);
end
z([1 end]) = zDomain;

% Sampled aligned family with the analytical barotropic sector.
h = solvedBasis.h(2:end);
FSampled = solvedBasis.F(z);
GSampled = solvedBasis.G(z);
HF = D*ones(1,nSolved+1);
HG = [0 D./h];
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
    [w,fitResidual] = fitBandWeights(FSampled(:,1:r+1),GSampled(:,2:r+1),HF(1:r+1),HG(2:r+1),GMetric,geometricWeights,D);
    row = bandDiagnostics(r,w,fitResidual,FSampled,GSampled,HF,HG,GMetric,options.gramTolerance,options.nOmittedModes,D, ...
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
        error("IMInternalModesDiscreteTransform:HydrostaticRootBandRejected", ...
            "The requested band r=%d fails: F energy error %.3g, G energy error %.3g (tolerance %.3g), F rank %d/%d, G rank %d/%d. Request fewer modes or omit nInternalModes for automatic selection.", ...
            last.internalModeCount,last.energyErrorF,last.energyErrorG,options.gramTolerance,last.rankF,last.internalModeCount+1,last.rankG,last.internalModeCount);
    end
    error("IMInternalModesDiscreteTransform:NoHydrostaticRootBand","No retained band passed the F and G energy tolerance %.3g on the %d-point root grid.",options.gramTolerance,nPoints);
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

gridDesign = struct(rule="hydrostaticRoot",generatingMode=nPoints-1,rootDesign=rootDesign, ...
    retainedInternalModeCount=r,weightRule="Nonnegative full-depth fit of stacked normalized F/G Gram errors for the retained band");
grid = IMVerticalGrid(z=z,weights=w,zDomain=zDomain,kind="modeRoots",design=gridDesign);
basis = solvedBasis.selectColumns(1:r+1);
transform = basis.discreteTransform(grid);
if options.forwardOperator == "truncatedCollocation"
    [forwardF,forwardG] = truncatedCollocation(r,collocationF,collocationG);
    transform = transform.transformWithForwardOperators(forwardF,forwardG,"truncatedCollocation");
end
assessment = IMHydrostaticRootTransformAssessment(selection=selection,nPoints=nPoints,retainedInternalModeCount=r, ...
    gramTolerance=options.gramTolerance,limitingCriterion=limitingCriterion,search=search,convergence=convergence, ...
    forwardOperator=options.forwardOperator,collocation=collocationDiagnostics, ...
    final=search(accepted,:));
end

function [basis,convergence,z,rootDesign] = convergedFamily(background,solver,refinedSolver,nSolved,nPoints,zDomain)
evp = IMInternalModes.hydrostaticGModes(background);
basis = IMHydrostaticModesBasis(solver.solveEVP(evp,nModes=nSolved),normalization="depthRMSF");
refined = IMHydrostaticModesBasis(refinedSolver.solveEVP(evp,nModes=nSolved),normalization="depthRMSF");
[zFine,wFine] = physicalReferenceQuadrature(zDomain);
labels = basis.columnLabels(2:end);
identity = struct(family="hydrostaticF",columnLabels=labels,normalization=basis.normalization,zDomain=zDomain);
candidate = basis.selectColumns(2:basis.nColumns);
reference = refined.selectColumns(2:refined.nColumns);
candidateModes = struct(identity=identity,values=struct(F=candidate.F(zFine),G=candidate.G(zFine)), ...
    equivalentDepths=candidate.h,provenance=struct(nEVP=solver.nEVP,coordinateKind=solver.coordinateKind));
referenceModes = struct(identity=identity,values=struct(F=reference.F(zFine),G=reference.G(zFine)), ...
    equivalentDepths=reference.h,provenance=struct(nEVP=refinedSolver.nEVP,coordinateKind=refinedSolver.coordinateKind));
measurements = assessModeConvergence(candidateModes,referenceModes,zFine,wFine).measurements;
[z,rootDesign] = basis.modeRootPoints(nPoints-1);
zRefined = refined.modeRootPoints(nPoints-1);
if numel(z) ~= nPoints || numel(zRefined) ~= nPoints
    error("IMInternalModesDiscreteTransform:HydrostaticRootCount","Generating G_%d must give %d points including endpoints; candidate/reference gave %d/%d.",nPoints-1,nPoints,numel(z),numel(zRefined));
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

function [w,residual] = fitBandWeights(F,G,HF,HG,GMetric,geometricWeights,depth)
% Stacked normalized Gram rows; off-diagonal rows carry sqrt(2) so the
% least-squares objective equals the full normalized Frobenius error.
[AF,bF] = gramRows(F,HF,ones(size(GMetric)));
[AG,bG] = gramRows(G,HG,GMetric);
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

function row = bandDiagnostics(r,w,fitResidual,F,G,HF,HG,GMetric,tolerance,nOmitted,depth,forwardOperator,collocationF,collocationG)
retainedF = 1:r+1;
retainedG = 2:r+1;
[epsF,rankF,condF,AF] = channelError(F(:,retainedF),w,HF(retainedF));
[epsG,rankG,condG,AG] = channelError(G(:,retainedG),GMetric.*w,HG(retainedG));
if forwardOperator == "truncatedCollocation"
    [AF,AGFull] = truncatedCollocation(r,collocationF,collocationG);
    AG = AGFull(2:end,:);
end
omitted = r+2:min(size(F,2),r+1+nOmitted);
if isempty(omitted)
    gainF = NaN; gainG = NaN;
else
    gainF = norm(sqrt(HF(retainedF)).'.*(AF*F(:,omitted))./sqrt(HF(omitted)),2);
    gainG = norm(sqrt(HG(retainedG)).'.*(AG*G(:,omitted))./sqrt(HG(omitted)),2);
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
    error("IMInternalModesDiscreteTransform:SingularHydrostaticRootCollocation","The square root-grid collocation system is numerically singular (reciprocal condition F %.3g, G %.3g).",rcondF,rcondG);
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
if isempty(H)
    epsilon = 0; rankValue = 0; conditionNumber = 1; forward = zeros(0,size(Phi,1));
    return
end
scale = 1./sqrt(H(:));
B = sqrt(metric).*Phi.*scale.';
Gamma = B.'*B;
Gamma = 0.5*(Gamma+Gamma.');
epsilon = norm(Gamma-eye(numel(H)),2);
singularValues = svd(B);
rankValue = sum(singularValues > max(size(B))*eps(max(singularValues)));
conditionNumber = cond(Gamma);
if rankValue < numel(H)
    forward = zeros(numel(H),size(Phi,1));
else
    [Q,R] = qr(B,0);
    forward = scale.*(R\(Q.'.*sqrt(metric).'));
end
end
