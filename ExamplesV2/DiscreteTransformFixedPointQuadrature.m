%% Scalar discrete transforms on fixed physical points
% This walkthrough constructs a Galerkin transform for hydrostatic G modes
% sampled on points chosen by the caller. The current discrete-transform API
% is scalar: because the EVP below uses the G formulation, the sampled scalar
% variable u and the transform basis Phi are G. Coupled F/G transforms are a
% later extension.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Solve a well-resolved exponential-stratification basis
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

nAvailableModes = 12;
nModes = 8;
nEVP = 160;
surfaceBoundary = IMBoundaryCondition.dirichlet();
bottomBoundary = IMBoundaryCondition.dirichlet();
evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g, ...
    surfaceBoundary=surfaceBoundary,bottomBoundary=bottomBoundary);
solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);
basisSet.normalization = "geostrophic";

%% Generate a deterministic mode-root candidate grid
% For the first `nModes` selected columns, `pointsFromModeRoots` returns both
% physical endpoints and the interior zeros of selected column `nModes+1`.
% The generating mode is already present here because the solve retained
% more than `nModes` columns. If it were absent, the basis set would solve
% for that one auxiliary mode automatically.
zModeRoot = basisSet.pointsFromModeRoots(nModes=nModes);
[~, modeRootFit] = basisSet.quadratureWeightsForPoints(z=zModeRoot,nModes=nModes);
fprintf("Mode-root grid: %d points; fitted Gram error: %.3e\n",length(zModeRoot),modeRootFit.transform.relativeGramError);

%% Supply fixed sample points and fit their weights
% These points include both boundaries and are refined toward the surface.
% `quadratureWeightsForPoints` retains the first `nModes` columns and chooses nonnegative
% weights that sum to the full depth.
nPoints = 24;
sigma = linspace(0,1,nPoints).';
z = zDomain(1) + D*(1 - (1 - sigma).^2);

[weights, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes);
transform = basisSet.discreteTransform(z=z,weights=weights,nModes=nModes);

% The ordinary one-line construction is equivalent:
% transform = basisSet.discreteTransform(z=z,nModes=nModes);
%
% Supplying weights bypasses fitting. For example, this reproduces the
% geometric comparison already stored in `fit`:
geometricTransform = basisSet.discreteTransform(z=z,weights=fit.geometricWeights,nModes=nModes);

%% Compare fitted and geometric quadrature
rule = ["fitted"; "geometric"];
relativeGramError = [transform.relativeGramError; geometricTransform.relativeGramError];
objectiveResidual = [fit.residualNorm; fit.geometricResidualNorm];
roundTripError = [transform.roundTripError; geometricTransform.roundTripError];
inverseMatrixConditionNumber = [transform.inverseMatrixConditionNumber; geometricTransform.inverseMatrixConditionNumber];
gramConditionNumber = [transform.gramConditionNumber; geometricTransform.gramConditionNumber];
depthSum = [sum(transform.weights); sum(geometricTransform.weights)];
minimumWeight = [min(transform.weights); min(geometricTransform.weights)];
maximumWeight = [max(transform.weights); max(geometricTransform.weights)];
diagnostics = table(rule,relativeGramError,objectiveResidual,roundTripError,inverseMatrixConditionNumber, ...
    gramConditionNumber,depthSum,minimumWeight,maximumWeight);

fprintf("\nFixed-point scalar transform for exponential hydrostatic G modes\n");
fprintf("Retained modes: %d; fixed points: %d; normalization: %s\n\n",nModes,nPoints,transform.normalization);
disp(diagnostics)

%% Transform known modal coefficients forward and back
% `transformBack` and `transformForward` accept multiple profile columns. A field made
% entirely from retained modes should round-trip to numerical precision.
mode = (1:nModes).';
coefficientsTrue = [1./mode (-1).^(mode - 1)./mode.^2];
valuesFromModes = transform.transformBack(coefficientsTrue);
coefficientsRecovered = transform.transformForward(valuesFromModes);
coefficientRoundTripError = norm(coefficientsRecovered - coefficientsTrue,2)/norm(coefficientsTrue,2);
fprintf("Relative coefficient round-trip error: %.3e\n",coefficientRoundTripError);

%% Transform a smooth sampled profile
% Galerkin projection finds the retained modal field whose residual is
% orthogonal to the retained basis in the sampled metric W.
profile = sin(pi*(z - zDomain(1))/D).*exp(z/900).*(1 + 0.2*cos(2*pi*(z - zDomain(1))/D));
profileCoefficients = transform.transformForward(profile);
profileReconstruction = transform.transformBack(profileCoefficients);
profileResidual = profile - profileReconstruction;
profileNorm = sqrt(profile.'*transform.metricMatrix*profile);
relativeProfileResidual = sqrt(profileResidual.'*transform.metricMatrix*profileResidual)/profileNorm;
fprintf("Relative sampled-metric residual of the smooth profile: %.3e\n",relativeProfileResidual);

%% Inspect the sampled modes, weights, and Gram errors
targetNorms = diag(transform.targetGramMatrix);
gramScale = 1./sqrt(abs(targetNorms));
fittedGramError = gramScale.*(transform.gramMatrix - transform.targetGramMatrix).*gramScale.';
geometricGramError = gramScale.*(geometricTransform.gramMatrix - geometricTransform.targetGramMatrix).*gramScale.';
gramColorLimit = max(abs([fittedGramError geometricGramError]),[],"all");

figure(Name="V2 fixed-point scalar discrete transform",Color="w");
tiledlayout(2,2,TileSpacing="compact",Padding="compact");

nexttile
plot(transform.inverseMatrix,z,LineWidth=1.1)
grid on
xlabel("G")
ylabel("z (m)")
title("Retained sampled modes")

nexttile
plot(fit.weights,z,"o-",LineWidth=1.1,MarkerSize=4)
hold on
plot(fit.geometricWeights,z,".-",LineWidth=1.1,MarkerSize=10)
hold off
grid on
xlabel("weight (m)")
ylabel("z (m)")
title("Quadrature weights")
legend(["fitted" "geometric"],Location="best")

nexttile
imagesc(transform.modeNumber,transform.modeNumber,geometricGramError)
axis xy image
clim([-gramColorLimit gramColorLimit])
colorbar
xlabel("mode number")
ylabel("mode number")
title("Geometric normalized Gram error")

nexttile
imagesc(transform.modeNumber,transform.modeNumber,fittedGramError)
axis xy image
clim([-gramColorLimit gramColorLimit])
colorbar
xlabel("mode number")
ylabel("mode number")
title("Fitted normalized Gram error")

figure(Name="V2 scalar discrete-transform profile projection",Color="w");
plot(profile,z,"k-",LineWidth=1.5)
hold on
plot(profileReconstruction,z,"--",LineWidth=1.5)
plot(profileResidual,z,":",LineWidth=1.2)
hold off
grid on
xlabel("amplitude")
ylabel("z (m)")
title("Smooth profile and retained-mode projection")
legend(["sampled profile" "reconstruction" "residual"],Location="best")
