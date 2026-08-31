%% Design a vertical grid for a hydrostatic Wave-Vortex basis
% Given a stratification profile and a fixed number of vertical points,
% this example determines where to place the points and how many vertical
% modes to retain. The design uses one grid and one set of integration
% weights for the paired F and G internal-mode structures.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Specify the stratification and point count
% The buoyancy frequency decreases exponentially with depth. We fix the
% number of physical points at Nz=64 and let the internal modes determine
% their locations.
D = 4000;
N0 = 5.2e-3;
b = 1300;
Nz = 128;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

%% Solve a hydrostatic F/G mode family
% The F formulation includes the barotropic mode. Its rigid-boundary
% conditions are equivalent to G=0 at the surface and bottom. Solving more
% modes than the point count gives the grid-design step enough modes to
% identify the complete candidate family. Both mode counts scale with Nz,
% so changing the requested point count does not impose a hidden limit.
nAvailableModes = Nz+8;
nEVP = max(256,4*Nz);
evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=zDomain);
solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=nAvailableModes);
basisSet.normalization = "geostrophic";

%% Choose the points and the retained modes
% `nPoints` is an exact physical point count. The returned grid contains
% both boundaries and the interior roots of the first mode beyond the
% candidate family. Quadrature weights are then fitted on those points for
% both F and G.
%
% Two error measures determine how many modes to keep:
%
% * Orthogonality error is the worst normalized difference between the
%   sampled and continuous mode inner products for F and G. The default
%   tolerance is 1e-2.
% * Quadratic-product error is the worst normalized projection error among
%   F_i F_j projected onto F, F_i G_j projected onto G, and G_i G_j
%   projected onto F. Here 0.1 is an illustrative nonlinear error budget.
quadraticProductTolerance = 0.1;
[transform,assessment] = basisSet.discreteTransform(nPoints=Nz,variables=["F","G"], ...
    quadraticAliasingTolerance=quadraticProductTolerance);

candidateModeCount = assessment.candidateModeCount;
linearModeCount = assessment.gramPolicy.maximumAcceptedModeCount;
quadraticModeCount = assessment.quadraticAliasingPolicy.maximumAcceptedModeCount;
retainedModeCount = assessment.retainedModeCount;
firstModeNumber = transform.modeNumber(1);
lastModeNumber = transform.modeNumber(end);

designSummary = table(Nz,candidateModeCount,linearModeCount,quadraticModeCount, ...
    retainedModeCount,firstModeNumber,lastModeNumber);
disp(designSummary)

%% Inspect where the points are placed
% The mode-root grid places more points where the vertical structures vary
% rapidly. For this profile that means fine spacing near the surface and
% progressively wider spacing at depth.
z = transform.z;
zIntegral = transform.weights;
zPlot = linspace(zDomain(1),zDomain(2),801).';
zMidpoint = 0.5*(z(1:end-1)+z(2:end));
dz = diff(z);

figure(Name="Wave-Vortex vertical grid",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(sqrt(N2(zPlot)),zPlot,"k-",LineWidth=1.5)
hold on
plot(sqrt(N2(z)),z,"o",MarkerSize=4,MarkerFaceColor="auto")
hold off
grid on
xlabel("N (s^{-1})")
ylabel("z (m)")
title(sprintf("%d mode-root points",Nz))
legend(["N(z)" "selected points"],Location="best")

nexttile
plot(dz,zMidpoint,"o-",LineWidth=1.2,MarkerSize=4)
grid on
xlabel("vertical spacing, \Delta z (m)")
ylabel("z (m)")
title("Grid spacing follows the stratification")

%% Inspect how the retained count is selected
% The curves use the same 64 points and fitted weights for every mode
% count. A count is accepted only while all smaller counts also satisfy the
% corresponding tolerance.
diagnostics = assessment.prefixDiagnostics;
modeCount = diagnostics.modeCount;
orthogonalityError = diagnostics.gramError;
quadraticProductError = diagnostics.quadraticAliasingError;

figure(Name="Wave-Vortex retained vertical modes",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
semilogy(modeCount,max(orthogonalityError,eps),"o-",LineWidth=1.2,MarkerSize=3,DisplayName="orthogonality error")
hold on
yline(assessment.gramPolicy.tolerance,"--",DisplayName="tolerance")
xline(linearModeCount,":",DisplayName=sprintf("retain through %d",linearModeCount))
hold off
grid on
xlabel("number of modes")
ylabel("normalized error")
title(sprintf("Linear transforms support %d modes",linearModeCount))
legend(Location="northwest")

nexttile
semilogy(modeCount,max(quadraticProductError,eps),"o-",LineWidth=1.2,MarkerSize=3,DisplayName="quadratic-product error")
hold on
yline(quadraticProductTolerance,"--",DisplayName="tolerance")
xline(quadraticModeCount,":",DisplayName=sprintf("retain through %d",quadraticModeCount))
hold off
ylim([1e-4 2])
grid on
xlabel("number of modes")
ylabel("normalized error")
title(sprintf("Nonlinear products support %d modes",quadraticModeCount))
legend(Location="northwest")

%% Extract the Wave-Vortex basis arrays
% These arrays contain the physical grid, integration weights, forward and
% inverse F/G transforms, physical mode numbers, and equivalent depths.
Finv = transform.inverseMatrix(variable="F");
F = transform.forwardMatrix(variable="F");
Ginv = transform.inverseMatrix(variable="G");
G = transform.forwardMatrix(variable="G");
j = transform.modeNumber(:);
h = transform.h(:);
Nj = length(j);

fprintf("The %d integration weights sum to %.6f m; minimum weight %.6f m.\n",length(zIntegral),sum(zIntegral),min(zIntegral));
fprintf("Retaining Nj=%d modes with physical mode numbers %d through %d.\n",Nj,j(1),j(end));

%% Verify a modal round trip
% G has no barotropic structure, so the first test coefficient is zero.
% The remaining coefficients can be reconstructed and recovered using
% either physical variable.
coefficient = exp(-(0:Nj-1).'/8).*cos((0:Nj-1).'*pi/5);
coefficient(1) = 0;
FRecovered = F*(Finv*coefficient);
GRecovered = G*(Ginv*coefficient);
relativeFError = norm(FRecovered-coefficient)/norm(coefficient);
relativeGError = norm(GRecovered(2:end)-coefficient(2:end))/norm(coefficient(2:end));

fprintf("F coefficient round-trip error: %.3e\n",relativeFError);
fprintf("Baroclinic G coefficient round-trip error: %.3e\n",relativeGError);
