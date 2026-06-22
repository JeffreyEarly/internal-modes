%% Surface geostrophic mode penetration depth in exponential stratification

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

N0 = 5.2e-3;
b = 1300;
f0 = 1e-4;
g0 = -N0*N0*b;
zDomain = [-4000 0];
N2 = @(z) N0*N0*exp(2*z/b);

k = logspace(log10(2*pi/1e6), log10(2*pi/10^2.7), 80);
L = 2*pi./k;
nEVP = 512;
z = linspace(zDomain(1), zDomain(2), 4000).';
threshold = 0.01;

problem = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=k, g0=g0);
solver = IMSolverSpectral(nEVP=nEVP);
basisSet = solver.solveSurfaceGeostrophicModes(problem);

h = basisSet.h;
iSignChange = find(h(1:end-1).*h(2:end) < 0);
kSignChange = zeros(size(iSignChange));
for iCross = 1:length(iSignChange)
    iLeft = iSignChange(iCross);
    fraction = -h(iLeft)/(h(iLeft+1) - h(iLeft));
    logKCrossing = log(k(iLeft)) + fraction*(log(k(iLeft+1)) - log(k(iLeft)));
    kSignChange(iCross) = exp(logKCrossing);
end

F = basisSet.F(z);
surfaceAmplitude = abs(F(end,:));
relativeAmplitude = abs(F)./surfaceAmplitude;
penetrationDepth = zeros(size(k));
fullDepth = zDomain(2) - zDomain(1);

for iMode = 1:length(k)
    iBelow = find(relativeAmplitude(:,iMode) <= threshold, 1, "last");
    if isempty(iBelow)
        penetrationDepth(iMode) = fullDepth;
        continue
    end

    if iBelow == length(z)
        penetrationDepth(iMode) = 0;
        continue
    end

    zBelow = z(iBelow);
    zAbove = z(iBelow+1);
    amplitudeBelow = relativeAmplitude(iBelow,iMode);
    amplitudeAbove = relativeAmplitude(iBelow+1,iMode);
    if amplitudeAbove == amplitudeBelow
        zCrossing = zBelow;
    else
        fraction = (threshold - amplitudeBelow)/(amplitudeAbove - amplitudeBelow);
        zCrossing = zBelow + fraction*(zAbove - zBelow);
    end
    penetrationDepth(iMode) = zDomain(2) - zCrossing;
end

figure(Name="V2 surface geostrophic mode penetration depth", Color="w");
tiledlayout(2, 1, TileSpacing="compact", Padding="compact");

topAxes = nexttile;
semilogx(k, h, LineWidth=1.5)
grid on
yline(0, "k-", LineWidth=0.75)
for iCross = 1:length(kSignChange)
    xline(kSignChange(iCross), "-", LineWidth=0.75, Color=[0.35 0.35 0.35])
end
ylabel("h (m)")
title("Surface-mode eigendepth and penetration depth")

bottomAxes = nexttile;
semilogx(k, -penetrationDepth, LineWidth=1.5)
grid on
for iCross = 1:length(kSignChange)
    xline(kSignChange(iCross), "-", LineWidth=0.75, Color=[0.35 0.35 0.35])
end
xlabel("\kappa (m^{-1})")
ylabel("z where |F(z)| < 0.01 |F(0)| (m)")
ylim(zDomain)
linkaxes([topAxes bottomAxes], "x")
