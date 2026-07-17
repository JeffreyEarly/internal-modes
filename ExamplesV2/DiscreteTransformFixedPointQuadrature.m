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

%% Supply fixed sample points and fit their increments
% These points include both boundaries and are refined toward the surface.
% `fitQuadrature` retains the first `nModes` columns and chooses nonnegative
% increments that sum to the full depth.
nPoints = 24;
sigma = linspace(0,1,nPoints).';
z = zDomain(1) + D*(1 - (1 - sigma).^2);

fit = basisSet.fitQuadrature(z=z,nModes=nModes);
transform = fit.fittedTransform;

% The ordinary one-line construction is equivalent:
% transform = basisSet.discreteTransform(z=z,nModes=nModes);
%
% Supplying increments bypasses fitting. For example, this reproduces the
% geometric comparison already stored in `fit`:
geometricTransform = basisSet.discreteTransform(z=z,increments=fit.geometricIncrements,nModes=nModes);

%% Compare fitted and geometric quadrature
rule = ["fitted"; "geometric"];
relativeGramError = [transform.relativeGramError; geometricTransform.relativeGramError];
objectiveResidual = [fit.fittedResidualNorm; fit.geometricResidualNorm];
roundTripError = [transform.roundTripError; geometricTransform.roundTripError];
basisConditionNumber = [transform.basisConditionNumber; geometricTransform.basisConditionNumber];
gramConditionNumber = [transform.gramConditionNumber; geometricTransform.gramConditionNumber];
depthSum = [sum(transform.increments); sum(geometricTransform.increments)];
minimumIncrement = [min(transform.increments); min(geometricTransform.increments)];
maximumIncrement = [max(transform.increments); max(geometricTransform.increments)];
diagnostics = table(rule,relativeGramError,objectiveResidual,roundTripError,basisConditionNumber, ...
    gramConditionNumber,depthSum,minimumIncrement,maximumIncrement);

fprintf("\nFixed-point scalar transform for exponential hydrostatic G modes\n");
fprintf("Retained modes: %d; fixed points: %d; normalization: %s\n\n",nModes,nPoints,transform.normalization);
disp(diagnostics)

%% Recover known modal coefficients
% `reconstruct` and `project` accept multiple profile columns. A field made
% entirely from retained modes should round-trip to numerical precision.
mode = (1:nModes).';
coefficientsTrue = [1./mode (-1).^(mode - 1)./mode.^2];
valuesFromModes = transform.reconstruct(coefficientsTrue);
coefficientsRecovered = transform.project(valuesFromModes);
coefficientRoundTripError = norm(coefficientsRecovered - coefficientsTrue,2)/norm(coefficientsTrue,2);
fprintf("Relative coefficient round-trip error: %.3e\n",coefficientRoundTripError);

%% Project a smooth sampled profile
% Galerkin projection finds the retained modal field whose residual is
% orthogonal to the retained basis in the sampled metric W.
profile = sin(pi*(z - zDomain(1))/D).*exp(z/900).*(1 + 0.2*cos(2*pi*(z - zDomain(1))/D));
profileCoefficients = transform.project(profile);
profileReconstruction = transform.reconstruct(profileCoefficients);
profileResidual = profile - profileReconstruction;
profileNorm = sqrt(profile.'*transform.metricMatrix*profile);
relativeProfileResidual = sqrt(profileResidual.'*transform.metricMatrix*profileResidual)/profileNorm;
fprintf("Relative sampled-metric residual of the smooth profile: %.3e\n",relativeProfileResidual);

%% Inspect the sampled modes, increments, and Gram errors
targetNorms = diag(transform.targetGramMatrix);
gramScale = 1./sqrt(abs(targetNorms));
fittedGramError = gramScale.*(transform.gramMatrix - transform.targetGramMatrix).*gramScale.';
geometricGramError = gramScale.*(geometricTransform.gramMatrix - geometricTransform.targetGramMatrix).*gramScale.';
gramColorLimit = max(abs([fittedGramError geometricGramError]),[],"all");

figure(Name="V2 fixed-point scalar discrete transform",Color="w");
tiledlayout(2,2,TileSpacing="compact",Padding="compact");

nexttile
plot(transform.basisMatrix,z,LineWidth=1.1)
grid on
xlabel("G")
ylabel("z (m)")
title("Retained sampled modes")

nexttile
plot(fit.fittedIncrements,z,"o-",LineWidth=1.1,MarkerSize=4)
hold on
plot(fit.geometricIncrements,z,".-",LineWidth=1.1,MarkerSize=10)
hold off
grid on
xlabel("increment (m)")
ylabel("z (m)")
title("Quadrature increments")
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
