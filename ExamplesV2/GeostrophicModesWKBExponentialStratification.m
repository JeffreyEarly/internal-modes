%% Positive-APV geostrophic modes with exponential stratification

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
g0 = -N0*N0*b;
zDomain = [-D 0];

N = @(z) N0*exp(z/b);
N2 = @(z) N0*N0*exp(2*z/b);

% D = 1300;
% N0 = 5.2e-3;
% g = 9.81;
% g0 = -N0*N0*D;
% zDomain = [-D 0];

N = @(z) N0*ones(size(z));
N2 = @(z) N0*N0*ones(size(z));

nModes = 4;
nEVP = 128;
z = linspace(zDomain(1), zDomain(2), 512).';

left = IMOperator() ...
    .plus(derivativeOrder=2) ...
    .plus(coefficient=@(z,~) -(2/b)*ones(size(z)), derivativeOrder=1);
right = IMOperator().plus(coefficient=@(z,ctx) -ctx.N2(z)/ctx.g, derivativeOrder=0);

surfaceLaw = IMOperator() ...
    .plus(derivativeOrder=1) ...
    .plus(coefficient=@(z,ctx) ctx.N2(z)*(1/ctx.g + 1/g0), derivativeOrder=0);
% surfaceBoundary = IMBoundary.custom(left=surfaceLaw, variable="F");
surfaceBoundary = IMBoundary.linearF(a=-(1/g + 1/g0), b=1);
bottomBoundary = IMBoundary.rigid();

innerWeights.F = @(z,~) ones(size(z));
innerWeights.G = @(z,ctx) ctx.N2(z)/ctx.g;
normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);

evp = IMEigenvalueProblem(name="wkbGeostrophicPositiveAPV", formulation="F", g=g, ...
    leftOperator=left, rightOperator=right, innerWeights=innerWeights, ...
    normalizations=normalizations, defaultNormalization=Normalization.unity, ...
    surfaceBoundary=surfaceBoundary, bottomBoundary=bottomBoundary, ...
    hFromEigenvalue=@(lambda) 1 ./ lambda, indexValidationMode="error");

solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
basisSet = solver.solveEVP(evp, nModes=nModes);

F = basisSet.F(z);
G = basisSet.G(z);
B = G - (1 + z/D).*F(end,:);

modeNumber = basisSet.modeNumber(:);
wkbIndex = modeNumber - 1;
integralN = N0*b*(1 - exp(-D/b));
theta0 = integralN ./ sqrt(g*basisSet.h(:));
anomalyRatio = max(abs(B), [], 1).' ./ abs(B(end,:)).';
leadingWkbRatio = leadingWkbAnomalyRatio(wkbIndex, D, b);

surfaceResidual = G(end,:) - (1 + g/g0)*F(end,:);
bottomResidual = G(1,:);
surfaceResidualScale = max(abs(G), [], 1);
bottomResidualScale = max(abs(G), [], 1);
surfaceResidualRelative = abs(surfaceResidual(:))./surfaceResidualScale(:);
bottomResidualRelative = abs(bottomResidual(:))./bottomResidualScale(:);

diagnostics = table(modeNumber, wkbIndex, basisSet.h(:), theta0, anomalyRatio, ...
    leadingWkbRatio, surfaceResidualRelative, bottomResidualRelative, ...
    VariableNames=["modeNumber" "wkbIndex" "h" "theta0" "anomalyRatio" ...
    "leadingWkbRatio" "surfaceBCRelative" "bottomBCRelative"]);

fprintf("Canonical exponential WKB estimate\n");
fprintf("D = %.0f m, N0 = %.3g s^-1, b = %.0f m, g0 = %.3g m s^-2, g/g0 = %.3g\n", ...
    D, N0, b, g0, g/g0);
fprintf("Lowest leading-WKB phase estimate without gradient corrections: theta0 = %.3f\n\n", ...
    sqrt((1 - exp(-D/b))*b*N0*N0*(1/g + 1/g0)));
disp(diagnostics)

modeLabels = "j = " + string(wkbIndex.');

figure(Name="V2 geostrophic positive-APV modes: exponential stratification", Color="w");
tiledlayout(1, 3, TileSpacing="compact", Padding="compact");

nexttile
plot(F, z, LineWidth=1.2)
grid on
xlabel("F")
ylabel("z (m)")
title("Solved F modes")
legend(modeLabels, Location="best")

nexttile
plot(G, z, LineWidth=1.2)
grid on
xlabel("G")
title("Diagnostic G modes")
legend(modeLabels, Location="best")

nexttile
plot(B, z, LineWidth=1.2)
grid on
xlabel("B")
title("Interior displacement anomaly")
legend(modeLabels, Location="best")

function ratio = leadingWkbAnomalyRatio(wkbIndex, D, b)
    ratio = ones(size(wkbIndex));
    baroclinicLike = wkbIndex >= 1;
    ratio(baroclinicLike) = wkbIndex(baroclinicLike) ...
        * pi*exp(D/(2*b))/(1 - exp(-D/b));
end
