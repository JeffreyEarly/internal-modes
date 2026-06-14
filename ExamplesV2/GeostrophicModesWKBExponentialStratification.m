%% Positive-APV geostrophic modes with exponential stratification

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
g0 = -N0*N0*b;
gd = -N0*N0*b/10;
zDomain = [-D 0];

N = @(z) N0*exp(z/b);
N2 = @(z) N0*N0*exp(2*z/b);

% D = 1300;
% N0 = 5.2e-3;
% g = 9.81;
% g0 = -N0*N0*D;
% zDomain = [-D 0];

% N = @(z) N0*ones(size(z));
% N2 = @(z) N0*N0*ones(size(z));

nModes = 4;
nEVP = 128;
z = linspace(zDomain(1), zDomain(2), 512).';

p = @(z,ctx) 1 ./ ctx.N2(z);
pz = @(z,ctx) -p(z,ctx) .* ctx.dzLogN2(z);

left = IMOperator() ...
    .plus(coefficient=p, derivativeOrder=2) ...
    .plus(coefficient=pz, derivativeOrder=1);

right = IMOperator() ...
    .plus(coefficient=@(z,ctx) -1/ctx.g, derivativeOrder=0);

% left = IMOperator() ...
%     .plus(derivativeOrder=2) ...
%     .plus(coefficient=@(z,~) -(2/b)*ones(size(z)), derivativeOrder=1);
% right = IMOperator().plus(coefficient=@(z,ctx) -ctx.N2(z)/ctx.g, derivativeOrder=0);

surfaceLaw = IMOperator() ...
    .plus(derivativeOrder=1) ...
    .plus(coefficient=@(z,ctx) ctx.N2(z)*(1/ctx.g + 1/g0), derivativeOrder=0);
% surfaceBoundary = IMBoundary.custom(left=surfaceLaw, variable="F");
surfaceBoundary = IMBoundary.linearF(a=-(1/g + 1/g0), b=1);
bottomBoundary = IMBoundary.linearF(a=1/gd, b=1);

innerWeights.F = @(z,~) ones(size(z));
innerWeights.G = @(z,ctx) ctx.N2(z)/ctx.g;
normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);

evp = IMEigenvalueProblem(name="unforced-APV-modes", formulation="F", g=g, ...
    leftOperator=left, rightOperator=right, innerWeights=innerWeights, ...
    normalizations=normalizations, defaultNormalization=Normalization.unity, ...
    surfaceBoundary=surfaceBoundary, bottomBoundary=bottomBoundary, ...
    hFromEigenvalue=@(lambda) 1 ./ lambda, indexValidationMode="error");

solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
basisSet = solver.solveEVP(evp, nModes=nModes);

F = basisSet.F(z);
G = basisSet.G(z);
B = G - (1 + z/D).*F(end,:);

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