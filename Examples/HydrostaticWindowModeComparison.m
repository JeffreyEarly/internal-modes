%% Compare hydrostatic window modes built from full F modes
% This example compares two ways to build partial-window modes from the full
% aligned hydrostatic F basis [F^0,F^1,...,F^N].
%
% The same aligned modal-coordinate vectors are then applied to the paired G
% basis [0,G^1,...,G^N]. The leading zero G column is the barotropic
% contribution: it is present in F and absent in G.
%
% The two window bases are:
%
% 1. visible modes ordered by h: find the visible subspace with the
%    full-depth-normalized F window operator, then diagonalize
%    diag(D,h_1,...,h_N) inside that visible subspace.
% 2. raw Gamma_win modes: diagonalize the unnormalized F window Gram matrix.
%    Its eigenvalues are partial-window kinetic-energy/equivalent-depth-like
%    norms in the named F basis.
% 3. generalized F/G modes: within the Boyd-retained raw subspace, solve a
%    generalized eigenvalue problem that diagonalizes both F-window and
%    G-window energy.

%% User settings
if ~exist('stratificationCase','var')
    stratificationCase = "exponential";
end
if ~exist('windowBounds','var')
    windowBounds = [-650 0];
end
if ~exist('nBaroclinicModes','var')
    nBaroclinicModes = 100;
end
if ~exist('nShow','var')
    nShow = 6;
end
if ~exist('nPlot','var')
    nPlot = 1024;
end
if ~exist('nQuad','var')
    nQuad = 4097;
end
if ~exist('visibilityThreshold','var')
    visibilityThreshold = 0.5;
end
if ~exist('minimumBaroclinicFractionForG','var')
    minimumBaroclinicFractionForG = 1e-3;
end
if ~exist('jStar','var')
    jStar = 3;
end
if ~exist('realizationSeed','var')
    realizationSeed = 42;
end

latitude = 31;
g = 9.81;

%% Build high-quality hydrostatic modes on plotting and window grids
modeData = buildHydrostaticModeData(stratificationCase,windowBounds,nBaroclinicModes,nPlot,nQuad,latitude,g);
nAlignedModes = nBaroclinicModes + 1;
nShow = min(nShow,nAlignedModes);

gammaF0 = [modeData.depth; modeData.h(:)];
sqrtGammaF0 = sqrt(gammaF0);
invSqrtGammaF0 = 1./sqrtGammaF0;

PhiFNormPlot = modeData.PhiFPlot .* invSqrtGammaF0.';
PhiFNormWindow = modeData.PhiFWindow .* invSqrtGammaF0.';
PhiGAlignedPlot = modeData.PhiGPlot;
PhiGAlignedWindow = modeData.PhiGWindow;

weightsFWindow = trapzWeights(modeData.zWindow);
weightsGWindow = weightsFWindow .* modeData.N2Window/modeData.g;

GammaFWin = modeData.PhiFWindow.'*(weightsFWindow.*modeData.PhiFWindow);
MWin = (invSqrtGammaF0.*GammaFWin).*invSqrtGammaF0.';
MWin = (MWin + MWin.')/2;
MGWin = PhiGAlignedWindow.'*(weightsGWindow.*PhiGAlignedWindow);
MGWin = (MGWin + MGWin.')/2;

%% Window basis 1: visible subspace ordered by h
[QVisibility,lambda] = sortedEigenvectors(MWin,"descend");
visibleIndex = lambda >= visibilityThreshold;
UVisible = QVisibility(:,visibleIndex);
if isempty(UVisible)
    error("HydrostaticWindowModeComparison:NoVisibleModes", ...
        "No window eigenvalues exceeded visibilityThreshold %.3g.", visibilityThreshold);
end
HAligned = diag(gammaF0);

HVisible = UVisible.'*HAligned*UVisible;
[SVisible,xi] = sortedEigenvectors(HVisible,"descend");
QH = orientColumns(UVisible*SVisible);
xi = xi(:);
lambdaEffH = diag(QH.'*MWin*QH);
barotropicFractionH = QH(1,:).^2;

%% Window basis 2: raw Gamma_win modes
[VRaw,muRaw] = sortedEigenvectors(GammaFWin,"descend");
VRaw = orientColumns(VRaw);
QRaw = sqrtGammaF0.*VRaw;
QRaw = QRaw ./ sqrt(sum(QRaw.^2,1));
QRaw = orientColumns(QRaw);

lambdaEffRaw = diag(QRaw.'*MWin*QRaw);
barotropicFractionRaw = QRaw(1,:).^2;
rawFullDepthGram = QRaw.'*QRaw;
rawWindowGram = VRaw.'*GammaFWin*VRaw;
[nProjectionRaw,nNegativeRaw,hasNegativeRaw] = boydRetainedCount(muRaw);
rawProjectionIndex = 1:nProjectionRaw;

%% Window basis 3: generalized F/G window modes
[QFG,lambdaFG,fgDiagnostics] = generalizedFGWindowBasis(MWin,MGWin,QRaw(:,rawProjectionIndex));
barotropicFractionFG = normalizedBarotropicFraction(QFG);

%% Build paired F and G modes from the same aligned eigenvectors
nShowH = min(nShow,size(QH,2));
nProjectionH = size(QH,2);
nShowRaw = min(nShow,length(muRaw));
rawShownIndex = 1:nShowRaw;
nProjectionFG = size(QFG,2);
nShowFG = min(nShow,nProjectionFG);

PsiFHPlot = PhiFNormPlot*QH(:,1:nShowH);
PsiFRawPlot = PhiFNormPlot*QRaw(:,rawShownIndex);
PsiFFGPlot = PhiFNormPlot*QFG(:,1:nShowFG);

PsiFHWindow = PhiFNormWindow*QH(:,1:nProjectionH);
PsiGHWindow = PhiGAlignedWindow*QH(:,1:nProjectionH);
PsiFRawWindow = PhiFNormWindow*QRaw(:,rawProjectionIndex);
PsiGRawWindow = PhiGAlignedWindow*QRaw(:,rawProjectionIndex);
PsiFFGWindow = PhiFNormWindow*QFG;
PsiGFGWindow = PhiGAlignedWindow*QFG;

PsiFHProjectionPlot = PhiFNormPlot*QH(:,1:nProjectionH);
PsiGHProjectionPlot = PhiGAlignedPlot*QH(:,1:nProjectionH);
PsiFRawProjectionPlot = PhiFNormPlot*QRaw(:,rawProjectionIndex);
PsiGRawProjectionPlot = PhiGAlignedPlot*QRaw(:,rawProjectionIndex);
PsiFFGProjectionPlot = PhiFNormPlot*QFG;
PsiGFGProjectionPlot = PhiGAlignedPlot*QFG;

gRetainedH = gActiveColumns(PsiGHWindow,weightsGWindow,QH(:,1:nProjectionH),minimumBaroclinicFractionForG);
gRetainedRaw = gActiveColumns(PsiGRawWindow,weightsGWindow,QRaw(:,rawProjectionIndex),minimumBaroclinicFractionForG);
gRetainedFG = gActiveColumns(PsiGFGWindow,weightsGWindow,QFG,minimumBaroclinicFractionForG);

PsiGHShownWindow = PhiGAlignedWindow*QH(:,1:nShowH);
PsiGRawShownWindow = PhiGAlignedWindow*QRaw(:,rawShownIndex);
PsiGFGShownWindow = PhiGAlignedWindow*QFG(:,1:nShowFG);
gShownH = gActiveColumns(PsiGHShownWindow,weightsGWindow,QH(:,1:nShowH),minimumBaroclinicFractionForG);
gShownRaw = gActiveColumns(PsiGRawShownWindow,weightsGWindow,QRaw(:,rawShownIndex),minimumBaroclinicFractionForG);
gShownFG = gActiveColumns(PsiGFGShownWindow,weightsGWindow,QFG(:,1:nShowFG),minimumBaroclinicFractionForG);
shownGIndexH = find(gShownH);
shownGIndexRaw = rawShownIndex(gShownRaw);
shownGIndexFG = find(gShownFG);
PsiGHPlot = PhiGAlignedPlot*QH(:,shownGIndexH);
PsiGRawPlot = PhiGAlignedPlot*QRaw(:,shownGIndexRaw);
PsiGFGPlot = PhiGAlignedPlot*QFG(:,shownGIndexFG);

%% Random realization with A/(j_*^2+j^2) spectrum
rng(realizationSeed,"twister");
targetSpectrum = targetModeSpectrum(nAlignedModes,jStar);
trueCoefficients = sqrt(targetSpectrum).*randn(nAlignedModes,1);

uPlot = PhiFNormPlot*trueCoefficients;
etaPlot = PhiGAlignedPlot*trueCoefficients;
uWindow = PhiFNormWindow*trueCoefficients;
etaWindow = PhiGAlignedWindow*trueCoefficients;

projectionH.F = projectWindowField(PsiFHWindow,weightsFWindow,uWindow);
projectionH.G = projectWindowField(PsiGHWindow(:,gRetainedH),weightsGWindow,etaWindow);
projectionRaw.F = projectWindowField(PsiFRawWindow,weightsFWindow,uWindow);
projectionRaw.G = projectWindowField(PsiGRawWindow(:,gRetainedRaw),weightsGWindow,etaWindow);
projectionFG.F = projectWindowField(PsiFFGWindow,weightsFWindow,uWindow);
projectionFG.G = projectWindowField(PsiGFGWindow(:,gRetainedFG),weightsGWindow,etaWindow);

uHPlot = PsiFHProjectionPlot*projectionH.F.coefficients;
etaHPlot = PsiGHProjectionPlot(:,gRetainedH)*projectionH.G.coefficients;
uRawPlot = PsiFRawProjectionPlot*projectionRaw.F.coefficients;
etaRawPlot = PsiGRawProjectionPlot(:,gRetainedRaw)*projectionRaw.G.coefficients;
uFGPlot = PsiFFGProjectionPlot*projectionFG.F.coefficients;
etaFGPlot = PsiGFGProjectionPlot(:,gRetainedFG)*projectionFG.G.coefficients;

projectionPlotIndex = true(size(modeData.zPlot));

%% Print diagnostics
fprintf('Hydrostatic window mode comparison\n');
fprintf('  stratification: %s\n', string(stratificationCase));
fprintf('  depth: %.0f m, window: [%.0f %.0f] m\n', modeData.depth, windowBounds(1), windowBounds(2));
fprintf('  aligned modes: %d = barotropic + %d baroclinic\n', nAlignedModes, nBaroclinicModes);
fprintf('  visible threshold %.3f retained %d visible directions\n', visibilityThreshold, size(QH,2));
fprintf('  G projection drops directions with baroclinic fraction <= %.1e\n', minimumBaroclinicFractionForG);
if hasNegativeRaw
    fprintf('  raw Gamma Boyd cutoff n_neg=%d retained %d directions\n', nNegativeRaw, nProjectionRaw);
else
    fprintf('  raw Gamma Boyd cutoff found no negative muRaw; retained %d of %d directions\n', ...
        nProjectionRaw, length(muRaw));
end
fprintf('  raw Gamma maximum visibility %.3f\n', max(lambdaEffRaw));
fprintf('  MWin eigenvalue bounds: %.3e\n', max([max(lambda - 1), max(-lambda), 0]));
fprintf('  visible-h basis orthogonality: %.3e\n', norm(QH.'*QH - eye(size(QH,2)),2));
fprintf('  raw Gamma_win diagonalization: %.3e\n', norm(rawWindowGram - diag(muRaw),2));
fprintf('  raw full-depth cross term after normalization: %.3e\n', offDiagonalNorm(rawFullDepthGram));
fprintf('  generalized F/G F diagonalization: %.3e\n', fgDiagnostics.FDiagonalizationError);
fprintf('  generalized F/G G diagonalization: %.3e\n', fgDiagnostics.GDiagonalizationError);
fprintf('  target spectrum A/(j_*^2+j^2), j_* = %.3g, sum = %.12f\n', jStar, sum(targetSpectrum));

fprintf('\nLeading visible-h modes\n');
printModeTable((1:nShowH).',xi(1:nShowH),lambdaEffH(1:nShowH),barotropicFractionH(1:nShowH), ...
    "xi","visibility");

fprintf('\nLeading raw Gamma_win modes\n');
printModeTable(rawShownIndex,muRaw(rawShownIndex),lambdaEffRaw(rawShownIndex), ...
    barotropicFractionRaw(rawShownIndex),"mu","visibility");

fprintf('\nLeading generalized F/G modes\n');
printModeTable((1:nShowFG).',lambdaFG(1:nShowFG),fgDiagnostics.FNorm(1:nShowFG), ...
    barotropicFractionFG(1:nShowFG),"G/F","F norm");

fprintf('\nProjection residual fractions\n');
fprintf('  visible-h F: %.3e, G: %.3e (using %d F modes, %d G-active modes; G rank %d/%d)\n', ...
    projectionH.F.residualFraction, projectionH.G.residualFraction, ...
    nProjectionH, nnz(gRetainedH), projectionH.G.rank, nnz(gRetainedH));
fprintf('  raw Gamma F: %.3e, G: %.3e (using %d F modes, %d G-active modes; G rank %d/%d)\n', ...
    projectionRaw.F.residualFraction, projectionRaw.G.residualFraction, ...
    nProjectionRaw, nnz(gRetainedRaw), projectionRaw.G.rank, nnz(gRetainedRaw));
fprintf('  generalized F/G F: %.3e, G: %.3e (using %d F modes, %d G-active modes; G rank %d/%d)\n', ...
    projectionFG.F.residualFraction, projectionFG.G.residualFraction, ...
    nProjectionFG, nnz(gRetainedFG), projectionFG.G.rank, nnz(gRetainedFG));

%% Figure 1: basis comparison
figure('Name',sprintf('Hydrostatic window modes: %s', string(stratificationCase)))
tiledlayout(3,3,'TileSpacing','compact','Padding','compact')

nexttile
stem(1:length(lambda),lambda,'filled','DisplayName','visibility \lambda'), hold on
stem(1:length(muRaw),muRaw/max(gammaF0),'DisplayName','raw \mu / max(\Gamma_0)')
yline(visibilityThreshold,'r--','visibility threshold','DisplayName','threshold')
xlabel('mode index')
ylabel('normalized value')
title('Window eigenvalues')
legend('Location','southwest')
ylim([-0.05 1.05])
grid on

nexttile
plot(1:length(lambdaEffH),lambdaEffH,'ko-','DisplayName','visible-h retained'), hold on
plot(1:length(lambdaEffRaw),lambdaEffRaw,'r.','DisplayName','raw candidates')
plot(rawProjectionIndex,lambdaEffRaw(rawProjectionIndex),'ro','DisplayName','raw Boyd retained')
yline(visibilityThreshold,'k--','threshold','DisplayName','threshold')
xlabel('basis mode index')
ylabel('window/full-depth F energy')
title('Effective visibility')
legend('Location','southwest')
ylim([-0.05 1.05])
grid on

nexttile
plotGeneralizedEigenvalues(lambdaFG)

nexttile
plotProfilesWithWindow(PsiFHPlot,modeData.zPlot,windowBounds, ...
    "F modes: visible subspace ordered by h",compose('%d: %.2g',(1:nShowH).',xi(1:nShowH)))

nexttile
plotProfilesWithWindow(PsiFRawPlot,modeData.zPlot,windowBounds, ...
    "F modes: raw \Gamma_{win} candidates",compose('%d: %.2g',rawShownIndex(:),muRaw(rawShownIndex)))

nexttile
plotProfilesWithWindow(PsiFFGPlot,modeData.zPlot,windowBounds, ...
    "F modes: generalized F/G",compose('%d: %.2g',(1:nShowFG).',lambdaFG(1:nShowFG)))

nexttile
plotProfilesWithWindow(PsiGHPlot,modeData.zPlot,windowBounds, ...
    "G modes: visible subspace ordered by h",compose('%d: %.2g',shownGIndexH(:),lambdaEffH(shownGIndexH)))

nexttile
plotProfilesWithWindow(PsiGRawPlot,modeData.zPlot,windowBounds, ...
    "G modes: raw \Gamma_{win} candidates",compose('%d: %.2g',shownGIndexRaw(:),lambdaEffRaw(shownGIndexRaw)))

nexttile
plotProfilesWithWindow(PsiGFGPlot,modeData.zPlot,windowBounds, ...
    "G modes: generalized F/G",compose('%d: %.2g',shownGIndexFG(:),lambdaFG(shownGIndexFG)))

%% Figure 2: modal coefficients
figure('Name',sprintf('Hydrostatic window coefficients: %s', string(stratificationCase)))
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

nexttile
plotCoefficientHeatmap(QH(:,1:nShowH),1:nShowH,0:nBaroclinicModes, ...
    "Visible-h coefficients","visible-h mode index")

nexttile
plotCoefficientHeatmap(QRaw(:,rawShownIndex),rawShownIndex,0:nBaroclinicModes, ...
    "Raw coefficients","raw \Gamma_{win} mode index")

nexttile
plotCoefficientHeatmap(QFG(:,1:nShowFG),1:nShowFG,0:nBaroclinicModes, ...
    "Generalized F/G coefficients","generalized mode index")

%% Figure 3: projection test
figure('Name',sprintf('Hydrostatic window projection: %s', string(stratificationCase)))
tiledlayout(3,2,'TileSpacing','compact','Padding','compact')

nexttile
plotProfileWithWindow(uPlot,modeData.zPlot,windowBounds,"random F realization")

nexttile
plotProfileWithWindow(etaPlot,modeData.zPlot,windowBounds,"paired G realization")

nexttile
plotWindowReconstruction(modeData.zPlot(projectionPlotIndex),uPlot(projectionPlotIndex), ...
    uHPlot(projectionPlotIndex),uRawPlot(projectionPlotIndex),uFGPlot(projectionPlotIndex), ...
    windowBounds,"F window reconstruction", ...
    [projectionH.F.residualFraction projectionRaw.F.residualFraction projectionFG.F.residualFraction])

nexttile
plotWindowReconstruction(modeData.zPlot(projectionPlotIndex),etaPlot(projectionPlotIndex), ...
    etaHPlot(projectionPlotIndex),etaRawPlot(projectionPlotIndex),etaFGPlot(projectionPlotIndex), ...
    windowBounds,"G window reconstruction", ...
    [projectionH.G.residualFraction projectionRaw.G.residualFraction projectionFG.G.residualFraction])

nexttile
stem(1:nProjectionH,projectionH.F.coefficients,'filled','DisplayName','visible-h F'), hold on
stem(rawProjectionIndex,projectionRaw.F.coefficients,'DisplayName','raw F')
stem(1:nProjectionFG,projectionFG.F.coefficients,'DisplayName','generalized F/G F')
xlabel('window mode index')
ylabel('coefficient')
title('Recovered F coefficients')
legend('Location','best')
grid on

nexttile
stem(find(gRetainedH),projectionH.G.coefficients,'filled','DisplayName','visible-h G'), hold on
stem(rawProjectionIndex(gRetainedRaw),projectionRaw.G.coefficients,'DisplayName','raw G')
stem(find(gRetainedFG),projectionFG.G.coefficients,'DisplayName','generalized F/G G')
xlabel('window mode index')
ylabel('coefficient')
title('Recovered G coefficients')
legend('Location','best')
grid on

drawnow

%% Local helpers
function data = buildHydrostaticModeData(stratificationCase,windowBounds,nBaroclinicModes,nPlot,nQuad,latitude,g)
caseName = string(stratificationCase);
switch caseName
    case "exponential"
        depth = 4000;
        N0 = 3*2*pi/3600;
        L_gm = 1300;
        N2Function = @(z) N0*N0*exp(2*z/L_gm);
        zDomain = [-depth 0];
        nEVP = max(256,ceil(2.5*(nBaroclinicModes + 1)));
        zPlot = linspace(zDomain(1),zDomain(2),nPlot).';
        zWindow = linspace(windowBounds(1),windowBounds(2),nQuad).';

        imPlot = InternalModesWKBSpectral(N2=N2Function,zIn=zDomain,zOut=zPlot, ...
            latitude=latitude,nEVP=nEVP,nModes=nBaroclinicModes,g=g);
        imPlot.normalization = Normalization.geostrophic;
        imPlot.upperBoundary = UpperBoundary.rigidLid;

        imWindow = InternalModesWKBSpectral(N2=N2Function,zIn=zDomain,zOut=zWindow, ...
            latitude=latitude,nEVP=nEVP,nModes=nBaroclinicModes,g=g);
        imWindow.normalization = Normalization.geostrophic;
        imWindow.upperBoundary = UpperBoundary.rigidLid;
    case "constant"
        depth = 1300;
        N0 = 5.2e-3;
        N2Function = @(z) N0*N0 + 0*z;
        zDomain = [-depth 0];
        zPlot = linspace(zDomain(1),zDomain(2),nPlot).';
        zWindow = linspace(windowBounds(1),windowBounds(2),nQuad).';

        imPlot = InternalModesConstantStratification(N0=N0,zIn=zDomain,zOut=zPlot, ...
            latitude=latitude,nModes=nBaroclinicModes,g=g);
        imPlot.upperBoundary = UpperBoundary.rigidLid;
        imPlot.normalization = Normalization.kConstant;

        imWindow = InternalModesConstantStratification(N0=N0,zIn=zDomain,zOut=zWindow, ...
            latitude=latitude,nModes=nBaroclinicModes,g=g);
        imWindow.upperBoundary = UpperBoundary.rigidLid;
        imWindow.normalization = Normalization.kConstant;
    otherwise
        error("HydrostaticWindowModeComparison:UnknownStratification", "Unknown stratificationCase ""%s"".", caseName);
end

if windowBounds(1) < zDomain(1) || windowBounds(2) > zDomain(2) || windowBounds(1) >= windowBounds(2)
    error("HydrostaticWindowModeComparison:InvalidWindow", ...
        "windowBounds must lie inside [%.3g %.3g].", zDomain(1), zDomain(2));
end

[FPlot,GPlot,h] = imPlot.modesAtFrequency(0);
[FWindow,GWindow] = imWindow.modesAtFrequency(0);

if caseName == "constant"
    geostrophicScale = sqrt((N0*N0 - imPlot.f0*imPlot.f0)/(N0*N0));
    FPlot = geostrophicScale*FPlot;
    GPlot = geostrophicScale*GPlot;
    FWindow = geostrophicScale*FWindow;
    GWindow = geostrophicScale*GWindow;
end

data.depth = depth;
data.g = g;
data.zPlot = zPlot;
data.zWindow = zWindow;
data.N2Plot = N2Function(zPlot);
data.N2Window = N2Function(zWindow);
data.h = h(1:nBaroclinicModes);
data.h = data.h(:);
data.PhiFPlot = cat(2,ones(size(zPlot)),FPlot(:,1:nBaroclinicModes));
data.PhiFWindow = cat(2,ones(size(zWindow)),FWindow(:,1:nBaroclinicModes));
data.PhiGPlot = cat(2,zeros(size(zPlot)),GPlot(:,1:nBaroclinicModes));
data.PhiGWindow = cat(2,zeros(size(zWindow)),GWindow(:,1:nBaroclinicModes));
end

function weights = trapzWeights(z)
weights = zeros(size(z));
weights(1) = (z(2) - z(1))/2;
weights(end) = (z(end) - z(end-1))/2;
weights(2:end-1) = (z(3:end) - z(1:end-2))/2;
end

function [Q,lambda] = sortedEigenvectors(A,sortDirection)
A = (A + A.')/2;
[Q,Lambda] = eig(A);
lambda = real(diag(Lambda));
[lambda,index] = sort(lambda,sortDirection);
Q = orientColumns(real(Q(:,index)));
end

function [Q,lambda,diagnostics] = generalizedFGWindowBasis(MF,MG,QBase)
if isempty(QBase)
    Q = zeros(size(MF,1),0);
    lambda = zeros(0,1);
    diagnostics.FDiagonalizationError = NaN;
    diagnostics.GDiagonalizationError = NaN;
    diagnostics.FNorm = zeros(0,1);
    return
end

MFSubspace = QBase.'*MF*QBase;
MGSubspace = QBase.'*MG*QBase;
MFSubspace = (MFSubspace + MFSubspace.')/2;
MGSubspace = (MGSubspace + MGSubspace.')/2;
[S,Lambda] = eig(MGSubspace,MFSubspace);
lambda = real(diag(Lambda));
[~,index] = sort(lambda,"descend");
Q = real(QBase*S(:,index));
Q = normalizeColumnsInMetric(Q,MF);
lambda = diag(Q.'*MG*Q);
Q = orientColumns(Q);

FGram = Q.'*MF*Q;
GGram = Q.'*MG*Q;
diagnostics.FDiagonalizationError = norm(FGram - eye(size(FGram)),2);
diagnostics.GDiagonalizationError = norm(GGram - diag(diag(GGram)),2);
diagnostics.FNorm = diag(FGram);
end

function Q = normalizeColumnsInMetric(Q,metric)
for iMode = 1:size(Q,2)
    normInMetric = sqrt(abs(Q(:,iMode).'*metric*Q(:,iMode)));
    Q(:,iMode) = Q(:,iMode)/normInMetric;
end
end

function [nRetained,nNegative,hasNegative] = boydRetainedCount(mu)
nNegative = find(mu < 0,1,'first');
hasNegative = ~isempty(nNegative);
if ~hasNegative
    nNegative = length(mu);
end
nRetained = floor(nNegative/2);
end

function fraction = normalizedBarotropicFraction(Q)
if isempty(Q)
    fraction = zeros(0,1);
    return
end
columnNorm = sum(Q.^2,1);
fraction = (Q(1,:).^2 ./ columnNorm).';
end

function isActive = gActiveColumns(PsiG,weightsG,modalCoordinates,minimumBaroclinicFraction)
if isempty(PsiG)
    isActive = false(1,0);
    return
end
columnNorms = sqrt(sum(weightsG.*PsiG.^2,1));
tolerance = max(size(PsiG))*eps(max(columnNorms));
coordinateNorm = sum(modalCoordinates.^2,1);
baroclinicFraction = sum(modalCoordinates(2:end,:).^2,1)./coordinateNorm;
isActive = columnNorms > tolerance & baroclinicFraction > minimumBaroclinicFraction;
end

function Q = orientColumns(Q)
for iMode = 1:size(Q,2)
    [~,iMax] = max(abs(Q(:,iMode)));
    Q(:,iMode) = sign(Q(iMax,iMode))*Q(:,iMode);
end
end

function spectrum = targetModeSpectrum(nModes,jStar)
modeNumbers = (0:(nModes-1)).';
shape = 1./(jStar*jStar + modeNumbers.^2);
spectrum = shape/sum(shape);
end

function result = projectWindowField(Psi,weights,field)
if isempty(Psi)
    result.coefficients = zeros(0,1);
    result.reconstruction = zeros(size(field));
    result.residualFraction = 1;
    result.rank = 0;
    result.conditionNumber = NaN;
    return
end
gram = Psi.'*(weights.*Psi);
rhs = Psi.'*(weights.*field);
tolerance = max(size(gram))*eps(norm(gram,2));
result.coefficients = pinv(gram,tolerance)*rhs;
result.reconstruction = Psi*result.coefficients;
residual = field - result.reconstruction;
result.residualFraction = weightedEnergy(residual,weights)/weightedEnergy(field,weights);
result.rank = rank(gram,tolerance);
result.conditionNumber = cond(gram);
end

function energy = weightedEnergy(field,weights)
energy = sum(weights.*field.^2);
end

function value = offDiagonalNorm(A)
value = norm(A - diag(diag(A)),'fro');
end

function printModeTable(modeIndex,values,metric,barotropicFraction,valueName,metricName)
fprintf('  %5s %12s %12s %12s\n', 'mode', valueName, metricName, 'q0^2');
for iMode = 1:length(values)
    fprintf('  %5d %12.4e %12.4e %12.4e\n', ...
        modeIndex(iMode), values(iMode), metric(iMode), barotropicFraction(iMode));
end
end

function plotGeneralizedEigenvalues(lambda)
if isempty(lambda)
    axis off
    text(0.5,0.5,"no generalized modes retained",HorizontalAlignment="center")
    title('Generalized F/G eigenvalues')
    return
end

stem(1:length(lambda),lambda,'filled')
xlabel('generalized mode index')
ylabel('G-window / F-window energy')
title('Generalized F/G eigenvalues')
grid on
end

function plotProfilesWithWindow(profiles,z,windowBounds,plotTitle,legendLabels)
if isempty(profiles)
    xLimits = [-1 1];
else
    xLimits = [min(profiles(:)) max(profiles(:))];
end
if diff(xLimits) == 0
    xLimits = xLimits + [-1 1];
end

patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
    [windowBounds(1) windowBounds(1) windowBounds(2) windowBounds(2)], ...
    [0.92 0.92 0.92],EdgeColor="none")
hold on
if isempty(profiles)
    text(0,mean(windowBounds),"no modes retained",HorizontalAlignment="center")
else
    plot(profiles,z,LineWidth=1.1)
    legend(legendLabels,Location="best")
end
xlabel('mode amplitude')
ylabel('z (m)')
title(plotTitle)
grid on
end

function plotProfileWithWindow(profile,z,windowBounds,plotTitle)
xLimits = [min(profile(:)) max(profile(:))];
if diff(xLimits) == 0
    xLimits = xLimits + [-1 1];
end

patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
    [windowBounds(1) windowBounds(1) windowBounds(2) windowBounds(2)], ...
    [0.92 0.92 0.92],EdgeColor="none")
hold on
plot(profile,z,'k',LineWidth=1.4)
xlabel('profile value')
ylabel('z (m)')
title(plotTitle)
grid on
end

function plotWindowReconstruction(z,trueProfile,visibleProfile,rawProfile,generalizedProfile, ...
    windowBounds,plotTitle,residuals)
xLimits = [min([trueProfile; visibleProfile; rawProfile; generalizedProfile]) ...
    max([trueProfile; visibleProfile; rawProfile; generalizedProfile])];
if diff(xLimits) == 0
    xLimits = xLimits + [-1 1];
end

patch([xLimits(1) xLimits(2) xLimits(2) xLimits(1)], ...
    [windowBounds(1) windowBounds(1) windowBounds(2) windowBounds(2)], ...
    [0.92 0.92 0.92],EdgeColor="none",DisplayName='window')
hold on
plot(trueProfile,z,'k',LineWidth=1.4,DisplayName='truth')
plot(visibleProfile,z,'b--',LineWidth=1.2,DisplayName='visible-h')
plot(rawProfile,z,'r:',LineWidth=1.4,DisplayName='raw \Gamma_{win}')
plot(generalizedProfile,z,'m-.',LineWidth=1.2,DisplayName='generalized F/G')
xlabel('profile value')
ylabel('z (m)')
title(sprintf('%s\nresiduals %.2e, %.2e, %.2e', plotTitle, residuals(1), residuals(2), residuals(3)))
legend(Location="best")
grid on
end

function plotCoefficientHeatmap(coefficients,xLabels,yLabels,plotTitle,xAxisLabel)
if isempty(coefficients)
    axis off
    text(0.5,0.5,"no modes retained",HorizontalAlignment="center")
    title(plotTitle)
    return
end

imagesc(1:length(xLabels),yLabels,coefficients)
set(gca,'YDir','normal')
set(gca,'XTick',1:length(xLabels),'XTickLabel',xLabels)
colorbar
xlabel(xAxisLabel)
ylabel('aligned modal index')
title(plotTitle)
end
