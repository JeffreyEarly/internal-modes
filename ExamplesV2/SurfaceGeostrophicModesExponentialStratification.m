%% Surface geostrophic modes in standard exponential stratification

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

N0 = 5.2e-3;
b = 1300;
f0 = 1e-4;
g0 = -N0*N0*b;
zDomain = [-4000 0];
N2 = @(z) N0*N0*exp(2*z/b);

L = [1e6 1e5 1e4 1e3];
k = 2*pi./L;
nEVP = 512;
z = linspace(zDomain(1), zDomain(2), 512).';

problem = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=k, g0=g0);
solver = IMSolverSpectral(nEVP=nEVP);
basisSet = solver.solveSurfaceGeostrophicModes(problem);

F = basisSet.F(z);
modeLabels = "L = " + string(L/1000) + " km, h = " + compose("%.3g", basisSet.h) + " m";

figure(Name="V2 surface geostrophic modes: exponential stratification", Color="w");
plot(F, z, LineWidth=1.2)
grid on
xlabel("F")
ylabel("z (m)")
title("Surface geostrophic modes")
legend(modeLabels, Location="best")
