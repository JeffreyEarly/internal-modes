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

N2 = @(z) N0*N0*exp(2*z/b);

nModes = 4;
nEVP = 128;
z = linspace(zDomain(1), zDomain(2), 512).';

p = @(z,ctx) 1 ./ ctx.N2(z);
q = @(z,~) zeros(size(z));
r = @(z,ctx) ones(size(z))/ctx.g;

surfaceBoundary = IMBoundaryCondition(a=-(1/g + 1/g0), b=1);
bottomBoundary = IMBoundaryCondition(a=1/gd, b=1);

normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);

evp = IMInternalModes(name="unforced-APV-modes", formulation="F", ...
    p=p, q=q, r=r, g=g, normalizations=normalizations, ...
    defaultNormalization=Normalization.unity, ...
    surfaceBoundary=surfaceBoundary, bottomBoundary=bottomBoundary);

solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, coordinateKind="wkb");
basisSet = solver.solveEVP(evp, nModes=nModes);

F = basisSet.F(z);
G = basisSet.G(z);
B = G - (1 + z/D).*F(end,:);

modeLabels = "j = " + string(basisSet.modeNumber);

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
