%% Custom quadrature objectives on fixed physical points
% The default fixed-point fit minimizes a normalized Gram residual. This
% walkthrough appends a small penalty on relative movement away from the
% geometric increments, illustrating how custom objectives can trade a small
% amount of Parseval accuracy for a less extreme quadrature rule.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Build the scalar hydrostatic G-mode basis
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

nAvailableModes = 12;
nModes = 8;
nPoints = 24;
evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g, ...
    surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
solver = IMSolverSpectral(nEVP=160,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);
basisSet.normalization = "geostrophic";

sigma = linspace(0,1,nPoints).';
z = zDomain(1) + D*(1 - (1 - sigma).^2);

%% Compute the default normalized-Gram fit
defaultFit = basisSet.fitQuadrature(z=z,nModes=nModes);

%% Add relative regularization toward geometric increments
% A custom callback receives the default normalized Gram system in
% `context.normalizedGramA` and `context.normalizedGramB`. Returning a new
% A and b replaces that objective, so callers may reweight existing rows or
% append new residuals. Here the additional rows penalize
%
%   sqrt(lambda) * (dz_i/dzGeometric_i - 1).
lambda = 1e-6;
objective = @(context) regularizedGramObjective(context,lambda);
regularizedFit = basisSet.fitQuadrature(z=z,nModes=nModes,objective=objective);

% Both fits retain the physical defaults `nonnegative=true` and
% `constrainDepth=true`. Algebraic investigations may disable them explicitly:
% algebraicFit = basisSet.fitQuadrature(z=z,nModes=nModes,nonnegative=false);
% freeDepthFit = basisSet.fitQuadrature(z=z,nModes=nModes,constrainDepth=false);

%% Compare every rule under the same Gram metric
% The custom objective contains extra regularization rows, so compare all
% candidates using the default normalized-Gram system rather than comparing
% their different full objective norms.
rule = ["geometric"; "normalized Gram"; "regularized Gram"];
increments = [defaultFit.geometricIncrements defaultFit.fittedIncrements regularizedFit.fittedIncrements];
gramResiduals = defaultFit.objectiveMatrix*increments - defaultFit.objectiveTarget;
normalizedGramResidual = vecnorm(gramResiduals,2,1).';
relativeGramError = [defaultFit.geometricTransform.relativeGramError; ...
    defaultFit.fittedTransform.relativeGramError; regularizedFit.fittedTransform.relativeGramError];
relativeIncrementDisplacement = vecnorm((increments - defaultFit.geometricIncrements)./defaultFit.geometricIncrements,2,1).'/sqrt(nPoints);
roundTripError = [defaultFit.geometricTransform.roundTripError; ...
    defaultFit.fittedTransform.roundTripError; regularizedFit.fittedTransform.roundTripError];
depthSum = sum(increments,1).';
minimumIncrement = min(increments,[],1).';
maximumIncrement = max(increments,[],1).';
diagnostics = table(rule,normalizedGramResidual,relativeGramError,relativeIncrementDisplacement, ...
    roundTripError,depthSum,minimumIncrement,maximumIncrement);

fprintf("\nCustom fixed-point quadrature objective\n");
fprintf("Relative geometric-increment regularization: lambda = %.1e\n\n",lambda);
disp(diagnostics)

%% Inspect the Parseval-versus-regularity tradeoff
plotLabels = ["geom" "Gram" "Gram + reg."];
figure(Name="V2 custom fixed-point quadrature objective",Color="w");
tiledlayout(1,3,TileSpacing="compact",Padding="compact");

nexttile
plot(increments(:,1),z,".-",LineWidth=1.1,MarkerSize=10)
hold on
plot(increments(:,2),z,"o-",LineWidth=1.1,MarkerSize=4)
plot(increments(:,3),z,"s-",LineWidth=1.1,MarkerSize=4)
hold off
grid on
xlabel("increment (m)")
ylabel("z (m)")
title("Increment rules")
legend(rule,Location="best")

nexttile
bar(relativeGramError)
set(gca,YScale="log",XTick=1:3,XTickLabel=plotLabels,XTickLabelRotation=20)
grid on
ylabel("relative Gram error")
title("Parseval accuracy")

nexttile
bar(max(relativeIncrementDisplacement,eps))
set(gca,YScale="log",XTick=1:3,XTickLabel=plotLabels,XTickLabelRotation=20)
grid on
ylabel("RMS relative displacement")
title("Distance from geometric increments")

%% Local objective builder
function specification = regularizedGramObjective(context,lambda)
relativeScale = 1./context.geometricIncrements;
regularizationMatrix = sqrt(lambda)*diag(relativeScale);
regularizationTarget = sqrt(lambda)*ones(length(context.z),1);
specification = struct("A",[context.normalizedGramA; regularizationMatrix], ...
    "b",[context.normalizedGramB; regularizationTarget], ...
    "name","normalizedGramWithGeometricRegularization");
end
