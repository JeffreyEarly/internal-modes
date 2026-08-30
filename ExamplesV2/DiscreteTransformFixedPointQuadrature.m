%% Fit a scalar discrete transform on caller-chosen physical points
% Suppose measurements already exist on a fixed vertical grid. This
% walkthrough fits quadrature weights on those points, builds a Galerkin
% transform for hydrostatic G modes, and separates the diagnostics that
% describe quadrature accuracy from those that describe linear algebra.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Solve exponential-stratification G modes
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

nAvailableModes = 12;
nModes = 8;
nEVP = 160;
evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g, ...
    surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);
basisSet.normalization = "geostrophic";

%% Choose a surface-refined grid and fit its quadrature weights
% The points include both boundaries and become more closely spaced toward
% the surface. Supplying `z` and `nModes` makes the requested eight-mode
% band strict: every requested mode must pass the enabled policies.
nPoints = 24;
sigma = linspace(0,1,nPoints).';
z = zDomain(1) + D*(1 - (1 - sigma).^2);

[transform, assessment] = basisSet.discreteTransform(z=z,nModes=nModes);
fit = assessment.weightFit;

% Supplying weights bypasses fitting. Use the geometric control-volume
% weights stored by the fit to build a comparison on exactly the same grid.
geometricTransform = basisSet.discreteTransform(z=z,weights=fit.geometricWeights,nModes=nModes);

%% Separate fit, transform, and conditioning diagnostics
% The fitted objective is the Frobenius norm of the normalized Gram error;
% it aggregates errors over all mode pairs. `relativeGramOperatorError` is
% the worst error over normalized modal combinations. `roundTripError`
% measures algebraic coefficient recovery and may remain tiny even when the
% quadrature rule has visible Gram error.
rule = ["fitted"; "geometric"];
frobeniusResidual = [fit.residualNorm; fit.geometricResidualNorm];
relativeGramOperatorError = [transform.relativeGramOperatorError; geometricTransform.relativeGramOperatorError];
roundTripError = [transform.roundTripError; geometricTransform.roundTripError];
forwardMatrixConditionNumber = [cond(transform.forwardMatrix); cond(geometricTransform.forwardMatrix)];
inverseMatrixConditionNumber = [transform.inverseMatrixConditionNumber; geometricTransform.inverseMatrixConditionNumber];
gramConditionNumber = [transform.gramConditionNumber; geometricTransform.gramConditionNumber];
depthSum = [sum(transform.weights); sum(geometricTransform.weights)];
minimumWeight = [min(transform.weights); min(geometricTransform.weights)];
maximumWeight = [max(transform.weights); max(geometricTransform.weights)];
diagnostics = table(rule,frobeniusResidual,relativeGramOperatorError,roundTripError, ...
    forwardMatrixConditionNumber,inverseMatrixConditionNumber,gramConditionNumber, ...
    depthSum,minimumWeight,maximumWeight);

fprintf("\nFixed-point transform for exponential hydrostatic G modes\n");
fprintf("Requested modes: %d; fixed points: %d; normalization: %s\n\n",nModes,nPoints,transform.normalization);
disp(diagnostics)

% Prefix diagnostics show how the same fitted rule behaves as modes are
% added. The Gram policy is the acceptance test; the other columns diagnose
% algebraic consistency and sensitivity.
prefixDiagnostics = assessment.prefixDiagnostics(:,["modeCount" "lastModeNumber" "gramError" ...
    "roundTripError" "inverseMatrixConditionNumber" "gramConditionNumber" "gramAccepted"]);
disp(prefixDiagnostics)

%% Recover known modal coefficients
% A field made entirely from retained modes should round-trip to numerical
% precision. Both transform directions accept multiple profile columns.
mode = (1:nModes).';
coefficientsTrue = [1./mode (-1).^(mode - 1)./mode.^2];
valuesFromModes = transform.transformBack(coefficientsTrue);
coefficientsRecovered = transform.transformForward(valuesFromModes);
coefficientRoundTripError = norm(coefficientsRecovered - coefficientsTrue,2)/norm(coefficientsTrue,2);
fprintf("Relative coefficient round-trip error: %.3e\n",coefficientRoundTripError);

%% Project a smooth sampled profile
% Galerkin projection returns the retained modal field whose residual is
% orthogonal to the retained basis in the sampled metric W. An arbitrary
% profile is therefore projected, not expected to round-trip exactly.
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

figure(Name="V2 fixed-point profile projection",Color="w");
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
