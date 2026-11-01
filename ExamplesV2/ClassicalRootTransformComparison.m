function results = ClassicalRootTransformComparison(options)
% Compare the classical root transform with the v4 algorithm on identical points.
%
% For each profile and point budget, `IMModes.classicalRootTransform` builds
% the root grid, fitted weights and qualified band. The v4 algorithm is
% reproduced from the same converged modes on the same points: square
% collocation with F = [1, F_1..F_{Nz-1}], the moment weights
% w4 = F^{-T}[D;0;...;0], and forward rows taken from the square inverse.
% Reported quantities:
%
% - usable band: the largest r with F and G energy errors <= tolerance under
%   the v4 moment weights and under the fitted weights (fitted per band);
% - omitted-mode gain: retained energy amplitude per unit energy amplitude
%   of internal modes r+1..r+8 under v4 collocation and the new Galerkin rows;
% - nonmodal input: Galerkin, prescribed-dual and continuous-projection F
%   coefficients of u(z) = exp(z/300), relative to the continuous projection;
% - timing: transform construction, and repeated batched forward/back
%   transforms of an Nz-by-nBatch array (dense multiplies in both cases).
%
% ```matlab
% results = ClassicalRootTransformComparison(nPoints=[17 33 65]);
% ```
%
% - Topic: Compare discrete transforms
% - Declaration: results = ClassicalRootTransformComparison(options)
% - Parameter options.nPoints: point budgets
% - Parameter options.gramTolerance: energy-error tolerance
% - Parameter options.nBatch: columns in the timing batch
% - Parameter options.outputFile: optional CSV path
% - Returns results: one table row per profile and point budget
arguments
    options.nPoints (1,:) double {mustBeInteger} = [17 33 65]
    options.gramTolerance (1,1) double = 1e-2
    options.nBatch (1,1) double {mustBeInteger,mustBePositive} = 128^2
    options.outputFile (1,1) string = ""
end
profiles = {
    "constant",       @(z) 2.5e-5*ones(size(z)),            1000
    "gmExponential",  @(z) (5.2e-3)^2*exp(2*z/1300),         4000
    "exponential400", @(z) 1e-4*exp(2*z/400),                1000
    "pycnocline",     @(z) 1e-6+4e-5*exp(-((z+300)/80).^2),  2000 };
rows = {};
for p = 1:size(profiles,1)
    [name,N2,D] = profiles{p,:};
    background = IMModeBackground(N2=N2,zDomain=[-D 0],latitude=33,surfaceCondition="rigidLid");
    for nPoints = options.nPoints
        t0 = tic;
        [transform,grid,assessment] = IMModes.classicalRootTransform(background,nPoints,gramTolerance=options.gramTolerance);
        setupTime = toc(t0);
        % Continuous reference modes at the converged resolution, beyond the retained band.
        solver = IMSolverSpectral(nEVP=2*assessment.convergence.candidateNEVP,coordinateKind="wkb");
        evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-D 0],latitude=33);
        nModes = nPoints+8;
        basis = solver.solveEVP(evp,nModes=nModes);
        z = grid.z; g = background.g; h = basis.h(1:nModes);
        F = [ones(nPoints,1) basis.F(z,normalization="geostrophic")];
        G = [zeros(nPoints,1) basis.G(z,normalization="geostrophic")]; G([1 end],:) = 0;
        HF = [D h]; GMetric = N2(z)/g;
        % v4 algorithm on the same points.
        Fsquare = F(:,1:nPoints);
        w4 = Fsquare.'\[D;zeros(nPoints-1,1)];
        collocationForward = inv(Fsquare); % v4 forms the explicit inverse
        usableV4 = 0;
        for r = nPoints-2:-1:1
            [eF,eG] = energyErrors(F(:,1:r+1),G(:,2:r+1),w4,HF(1:r+1),GMetric);
            if eF <= options.gramTolerance && eG <= options.gramTolerance, usableV4 = r; break, end
        end
        r = assessment.retainedInternalModeCount;
        [eF4,eG4] = energyErrors(F(:,1:r+1),G(:,2:r+1),w4,HF(1:r+1),GMetric);
        omitted = r+2:min(r+9,nModes+1);
        galerkinF = transform.forwardMatrix(variable="F");
        gainNew = norm(sqrt(HF(1:r+1)).'.*(galerkinF*F(:,omitted))./sqrt(HF(omitted)),2);
        v4Rows = collocationForward(1:r+1,:);
        gainV4 = norm(sqrt(HF(1:r+1)).'.*(v4Rows*F(:,omitted))./sqrt(HF(omitted)),2);
        % Nonmodal input u = exp(z/300).
        fine = IMVerticalGrid.physicalDepth([-D 0],2049);
        Ffine = [ones(fine.nPoints,1) basis.F(fine.z,normalization="geostrophic")];
        u = @(zz) exp(zz/300);
        continuous = (Ffine(:,1:r+1).'*(fine.weights.*u(fine.z)))./HF(1:r+1).';
        galerkin = galerkinF*u(z);
        dual = transform.modeProjectionFunctional(u(z),variable="F")./HF(1:r+1).';
        collocation = v4Rows*u(z);
        relative = @(c) norm(sqrt(HF(1:r+1)).'.*(c-continuous))/norm(sqrt(HF(1:r+1)).'.*continuous);
        % Timing of repeated batched transforms.
        X = randn(nPoints,options.nBatch);
        GForward = transform.forwardMatrix(variable="G"); FBack = transform.inverseMatrix(variable="F");
        forwardTime = timeit(@() galerkinF*X);
        backTime = timeit(@() FBack*(GForward*X));
        v4Time = timeit(@() v4Rows*X);
        rows(end+1,:) = {name,nPoints,r,usableV4,assessment.final.energyErrorF,assessment.final.energyErrorG,eF4,eG4, ...
            gainNew,gainV4,relative(galerkin),relative(dual),relative(collocation),min(grid.weights),min(w4), ...
            assessment.convergence.candidateNEVP,setupTime,forwardTime,backTime,v4Time}; %#ok<AGROW>
    end
end
results = cell2table(rows,VariableNames=["profile","nPoints","retainedNew","usableV4Weights","energyErrorFNew","energyErrorGNew", ...
    "energyErrorFV4WeightsSameBand","energyErrorGV4WeightsSameBand","omittedGainFNew","omittedGainFV4Collocation", ...
    "nonmodalErrorGalerkin","nonmodalErrorPrescribedDual","nonmodalErrorV4Collocation","minimumWeightNew","minimumWeightV4", ...
    "nEVP","setupSeconds","forwardSeconds","forwardBackSeconds","v4ForwardSeconds"]);
if options.outputFile ~= ""
    writetable(results,options.outputFile);
end
end

function [epsF,epsG] = energyErrors(F,G,w,HF,GMetric)
gram = F.'*(w.*F);
epsF = norm((gram-diag(HF))./sqrt(HF(:)*HF(:).'),2);
gram = G.'*((GMetric.*w).*G);
epsG = norm(gram-eye(size(G,2)),2);
end
