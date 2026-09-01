%% Design a vertical grid for free-surface Wave-Vortex modes
% Given a stratification profile and a fixed number of vertical points,
% this example builds geostrophic available-potential-vorticity (APV) and
% mean-density-anomaly (MDA) transforms on one physical grid. It compares
% the APV F/G Gram errors, the MDA G Gram error, and the APV coupled
% quadratic-product error that limit their independently retained bands.

repoRoot = fileparts(fileparts(mfilename("fullpath")));
addpath(repoRoot);

%% Specify the stratification, endpoints, and point count
% The buoyancy frequency decreases exponentially with depth. We fix the
% number of physical points in Nz and let the APV modes determine their
% locations. Both families then fit their own quadrature weights on those
% shared points. Positive finite endpoint accelerations give both families
% positive-definite projection metrics and keep both endpoints active.
D = 4000;
N0 = 5.2e-3;
b = 1300;
g = 9.81;
g0 = 0.02;
gd = 0.03;
Nz = 128;
zDomain = [-D 0];
N2 = @(z) N0*N0*exp(2*z/b);

%% Solve the APV and MDA mode families
% Solving slightly more modes than the point count gives the APV grid
% designer enough modes to identify its complete candidate family. The MDA
% solve uses the same numerical resolution but remains a separate modal
% family with its own retained count.
nAvailableModes = Nz+8;
nEVP = max(256,4*Nz);
solver = IMSolverSpectral(nEVP=nEVP,coordinateKind="wkb");

apvEVP = IMInternalModes.geostrophicAPVModes( ...
    N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd,surfaceBoundary="freeSurface");
apvBasis = solver.solveEVP(apvEVP,nModes=nAvailableModes);

mdaEVP = IMInternalModes.meanDensityAnomalyModes( ...
    N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd);
mdaBasis = solver.solveEVP(mdaEVP,nModes=nAvailableModes);

%% Build both transforms on the APV-designed grid
% `nPoints` is an exact physical point count. The APV transform first
% designs a mode-root grid, fits weights, and assesses every candidate
% prefix. The MDA transform then receives those points explicitly and fits
% an independent G-channel quadrature rule. The families are not forced to
% retain the same number of modes.
%
% The Gram error measures normalized distortion of the continuous modal
% inner product by the sampled transform. APV reports both F and G channel
% errors and applies their worst value. MDA is projected through G, while
% its aligned F structures remain available for pressure synthesis.
%
% APV additionally assesses the worst projection error in the three
% coupled product channels
%
%   F_i F_j -> F,   F_i G_j -> G,   G_i G_j -> F.
%
% MDA does not use this APV quadratic-product policy.
gramTolerance = 1e-2;
quadraticAliasingTolerance = 0.1;
[apvTransform,apvAssessment] = apvBasis.discreteTransform( ...
    nPoints=Nz,variables=["F","G"],gramTolerance=gramTolerance, ...
    quadraticAliasingTolerance=quadraticAliasingTolerance);

z = apvTransform.z;
[mdaTransform,mdaAssessment] = mdaBasis.discreteTransform( ...
    z=z,variables="G",gramTolerance=gramTolerance);
[mdaDesignedTransform,mdaDesignedAssessment] = mdaBasis.discreteTransform( ...
    nPoints=Nz,variables="G",gramTolerance=gramTolerance);

%% Summarize the independently retained modal bands
% `candidateModeCount` is the complete band tested on the fitted rule.
% `maximumAcceptedModeCount` is the largest cumulative prefix satisfying a
% policy. The production transform retains the intersection of the enabled
% policies for that family. The third row is a diagnostic comparison: it
% lets MDA design its own point locations instead of using the APV grid.
family = ["APV";"MDA";"MDA"];
pointRule = ["APV-designed";"APV-designed";"MDA-designed"];
candidateModeCount = [apvAssessment.candidateModeCount; ...
    mdaAssessment.candidateModeCount;mdaDesignedAssessment.candidateModeCount];
gramAcceptedModeCount = [apvAssessment.gramPolicy.maximumAcceptedModeCount; ...
    mdaAssessment.gramPolicy.maximumAcceptedModeCount; ...
    mdaDesignedAssessment.gramPolicy.maximumAcceptedModeCount];
quadraticAcceptedModeCount = [apvAssessment.quadraticAliasingPolicy.maximumAcceptedModeCount;NaN;NaN];
retainedModeCount = [apvAssessment.retainedModeCount; ...
    mdaAssessment.retainedModeCount;mdaDesignedAssessment.retainedModeCount];
retainedNegativeModeCount = [nnz(apvTransform.modeNumber < 0); ...
    nnz(mdaTransform.modeNumber < 0);nnz(mdaDesignedTransform.modeNumber < 0)];
retainsZeroMode = [any(apvTransform.modeNumber == 0); ...
    any(mdaTransform.modeNumber == 0);any(mdaDesignedTransform.modeNumber == 0)];
firstRetainedModeLabel = [apvTransform.modeNumber(1); ...
    mdaTransform.modeNumber(1);mdaDesignedTransform.modeNumber(1)];
lastRetainedModeLabel = [apvTransform.modeNumber(end); ...
    mdaTransform.modeNumber(end);mdaDesignedTransform.modeNumber(end)];

designSummary = table(family,pointRule,candidateModeCount,gramAcceptedModeCount, ...
    quadraticAcceptedModeCount,retainedModeCount,retainedNegativeModeCount, ...
    retainsZeroMode,firstRetainedModeLabel,lastRetainedModeLabel);
disp(designSummary)

apvDiagnostics = apvAssessment.prefixDiagnostics;
mdaDiagnostics = mdaAssessment.prefixDiagnostics;
mdaDesignedDiagnostics = mdaDesignedAssessment.prefixDiagnostics;
firstAPVQuadraticRejection = find(~apvDiagnostics.quadraticAccepted,1);
if ~isempty(firstAPVQuadraticRejection)
    rejected = apvDiagnostics(firstAPVQuadraticRejection,:);
    fprintf("First rejected APV quadratic prefix: %d modes; %s from physical modes %g and %g; error %.3e.\n", ...
        rejected.modeCount,rejected.quadraticLimitingChannel, ...
        rejected.quadraticLimitingModeNumberI,rejected.quadraticLimitingModeNumberJ, ...
        rejected.quadraticAliasingError);
end

firstMDAGramRejection = find(~mdaDiagnostics.gramAccepted,1);
if ~isempty(firstMDAGramRejection)
    rejected = mdaDiagnostics(firstMDAGramRejection,:);
    fprintf("First rejected MDA Gram prefix: %d modes; limiting channel %s; error %.3e.\n", ...
        rejected.modeCount,rejected.gramLimitingVariable,rejected.gramError);
end

firstMDADesignedGramRejection = find(~mdaDesignedDiagnostics.gramAccepted,1);
if ~isempty(firstMDADesignedGramRejection)
    rejected = mdaDesignedDiagnostics(firstMDADesignedGramRejection,:);
    fprintf("First MDA-designed Gram rejection: %d modes; limiting channel %s; error %.3e.\n", ...
        rejected.modeCount,rejected.gramLimitingVariable,rejected.gramError);
end

%% Compare the APV- and MDA-designed point rules
% The APV mode-root grid places more points where the vertical structures
% vary rapidly. For this profile that means fine spacing near the surface
% and progressively wider spacing at depth. The second rule shows where
% MDA would place the same number of points if it designed its own grid.
zPlot = linspace(zDomain(1),zDomain(2),801).';
mdaDesignedZ = mdaDesignedTransform.z;
apvZMidpoint = 0.5*(z(1:end-1)+z(2:end));
mdaZMidpoint = 0.5*(mdaDesignedZ(1:end-1)+mdaDesignedZ(2:end));
apvDz = diff(z);
mdaDz = diff(mdaDesignedZ);

figure(Name="Free-surface Wave-Vortex vertical grid",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
plot(sqrt(N2(zPlot)),zPlot,"k-",LineWidth=1.5)
hold on
plot(sqrt(N2(z)),z,"o",MarkerSize=4,DisplayName="APV-designed points")
plot(sqrt(N2(mdaDesignedZ)),mdaDesignedZ,"x",MarkerSize=5, ...
    LineWidth=1.0,DisplayName="MDA-designed points")
hold off
grid on
xlabel("N (s^{-1})")
ylabel("z (m)")
title(sprintf("Two family-designed %d-point grids",Nz))
legend(["N(z)" "APV-designed points" "MDA-designed points"],Location="best")

nexttile
plot(apvDz,apvZMidpoint,"o-",LineWidth=1.2,MarkerSize=4,DisplayName="APV-designed")
hold on
plot(mdaDz,mdaZMidpoint,"x-",LineWidth=1.2,MarkerSize=4,DisplayName="MDA-designed")
hold off
grid on
xlabel("vertical spacing, \Delta z (m)")
ylabel("z (m)")
title("Family-designed grid spacing")
legend(Location="best")

%% Inspect the lowest modes on both point rules
% The lines show each G mode on a fine vertical grid. Circles mark the
% APV-designed points used by the shared-grid transform, while crosses in
% the MDA panel mark the independently MDA-designed points. Each mode is
% divided by its own maximum absolute value so that the sampling of every
% shape remains visible. The basis already applies the canonical
% shallow-interior-G-positive orientation, so the example does not apply a
% separate display phase.
%
% Array columns are ordinal storage locations, while `modeNumber` is the
% physical label. Negative labels identify retained negative eigenvalues,
% zero identifies an actual null mode, and positive labels identify the
% ordinary positive branch. APV and MDA labels are separate coordinates;
% the family names are therefore included explicitly in the legends. The
% first four MDA candidates are shown even when the APV-designed rule
% accepts only a shorter MDA prefix.
nModeShapes = 4;
apvModeShapeCount = min(nModeShapes,length(apvBasis.modeNumber));
mdaModeShapeCount = min(nModeShapes,length(mdaBasis.modeNumber));

apvGModeFine = apvBasis.G(zPlot);
apvGModeFine = apvGModeFine(:,1:apvModeShapeCount);
apvGModeScale = max(abs(apvGModeFine),[],1);
apvGModeFine = apvGModeFine./apvGModeScale;
apvGModeAtPoints = apvBasis.G(z);
apvGModeAtPoints = apvGModeAtPoints(:,1:apvModeShapeCount)./apvGModeScale;

mdaGModeFine = mdaBasis.G(zPlot);
mdaGModeFine = mdaGModeFine(:,1:mdaModeShapeCount);
mdaGModeScale = max(abs(mdaGModeFine),[],1);
mdaGModeFine = mdaGModeFine./mdaGModeScale;
mdaGModeAtPoints = mdaBasis.G(z);
mdaGModeAtPoints = mdaGModeAtPoints(:,1:mdaModeShapeCount)./mdaGModeScale;
mdaGModeAtDesignedPoints = mdaBasis.G(mdaDesignedZ);
mdaGModeAtDesignedPoints = mdaGModeAtDesignedPoints(:,1:mdaModeShapeCount)./mdaGModeScale;

modeColors = colororder;
zeroModeColor = [0.35 0.35 0.35];
negativeModeColors = [0.64 0.08 0.18;0.85 0.33 0.10];
figure(Name="Lowest free-surface Wave-Vortex modes",Color="w");
tiledlayout(1,2,TileSpacing="compact",Padding="compact");

nexttile
hold on
for iMode = 1:apvModeShapeCount
    modeNumber = apvBasis.modeNumber(iMode);
    if modeNumber < 0
        modeKind = "negative";
        modeColor = negativeModeColors(min(abs(modeNumber),size(negativeModeColors,1)),:);
    elseif modeNumber == 0
        modeKind = "zero/null";
        modeColor = zeroModeColor;
    else
        modeKind = "positive";
        colorIndex = mod(modeNumber-1,size(modeColors,1))+1;
        modeColor = modeColors(colorIndex,:);
    end
    plot(apvGModeFine(:,iMode),zPlot,LineWidth=1.4, ...
        Color=modeColor,DisplayName="apvModeNumber="+string(modeNumber)+" ("+modeKind+")")
    plot(apvGModeAtPoints(:,iMode),z,"o",Color=modeColor, ...
        MarkerSize=3,MarkerFaceColor="w",HandleVisibility="off")
end
plot(NaN,NaN,"ko",MarkerSize=4,MarkerFaceColor="w",DisplayName="APV-designed points")
hold off
grid on
xlabel("canonically oriented, individually normalized G_j(z)")
ylabel("z (m)")
title("Lowest APV modes")
legend(Location="best")

nexttile
hold on
for iMode = 1:mdaModeShapeCount
    modeNumber = mdaBasis.modeNumber(iMode);
    if modeNumber < 0
        modeKind = "negative";
        modeColor = negativeModeColors(min(abs(modeNumber),size(negativeModeColors,1)),:);
    elseif modeNumber == 0
        modeKind = "zero/null";
        modeColor = zeroModeColor;
    else
        modeKind = "positive";
        colorIndex = mod(modeNumber-1,size(modeColors,1))+1;
        modeColor = modeColors(colorIndex,:);
    end
    plot(mdaGModeFine(:,iMode),zPlot,LineWidth=1.4, ...
        Color=modeColor,DisplayName="mdaModeNumber="+string(modeNumber)+" ("+modeKind+")")
    plot(mdaGModeAtPoints(:,iMode),z,"o",Color=modeColor, ...
        MarkerSize=3,MarkerFaceColor="w",HandleVisibility="off")
    plot(mdaGModeAtDesignedPoints(:,iMode),mdaDesignedZ,"x",Color=modeColor, ...
        MarkerSize=4,LineWidth=1.0,HandleVisibility="off")
end
plot(NaN,NaN,"ko",MarkerSize=4,MarkerFaceColor="w",DisplayName="APV-designed points")
plot(NaN,NaN,"kx",MarkerSize=5,LineWidth=1.0,DisplayName="MDA-designed points")
hold off
grid on
xlabel("canonically oriented, individually normalized G_j(z)")
ylabel("z (m)")
title(["Lowest MDA modes",sprintf("%d pass APV points; %d pass MDA points", ...
    mdaAssessment.gramPolicy.maximumAcceptedModeCount, ...
    mdaDesignedAssessment.gramPolicy.maximumAcceptedModeCount)])
legend(Location="best")

%% Compare the Gram and APV quadratic errors
% All curves assess leading modal prefixes without refitting their family's
% weights. The top-left panel separates APV F and G errors so the limiting
% physical variable is visible. The top-right panel shows why the same
% points can support a different MDA band. The lower panels show the APV
% nonlinear limit and the two independently fitted weight vectors.
apvFDiagnostics = apvAssessment.variablePrefixDiagnostics(variable="F");
apvGDiagnostics = apvAssessment.variablePrefixDiagnostics(variable="G");
mdaGDiagnostics = mdaAssessment.variablePrefixDiagnostics(variable="G");
mdaDesignedGDiagnostics = mdaDesignedAssessment.variablePrefixDiagnostics(variable="G");

figure(Name="Free-surface Wave-Vortex transform diagnostics",Color="w");
tiledlayout(2,2,TileSpacing="compact",Padding="compact");

nexttile
semilogy(apvFDiagnostics.modeCount,max(apvFDiagnostics.gramError,eps),"o-", ...
    LineWidth=1.1,MarkerSize=3,DisplayName="F Gram error")
hold on
semilogy(apvGDiagnostics.modeCount,max(apvGDiagnostics.gramError,eps),"o-", ...
    LineWidth=1.1,MarkerSize=3,DisplayName="G Gram error")
semilogy(apvDiagnostics.modeCount,max(apvDiagnostics.gramError,eps),"k-", ...
    LineWidth=1.5,DisplayName="worst channel")
yline(gramTolerance,"--",DisplayName="tolerance")
xline(apvAssessment.gramPolicy.maximumAcceptedModeCount,":", ...
    DisplayName=sprintf("accept through %d",apvAssessment.gramPolicy.maximumAcceptedModeCount))
hold off
grid on
xlabel("APV prefix mode count")
ylabel("normalized Gram error")
title("APV F/G Gram errors")
legend(Location="northwest")

nexttile
semilogy(mdaGDiagnostics.modeCount,max(mdaGDiagnostics.gramError,eps),"o-", ...
    LineWidth=1.1,MarkerSize=3,DisplayName="APV-designed points")
hold on
semilogy(mdaDesignedGDiagnostics.modeCount,max(mdaDesignedGDiagnostics.gramError,eps),"x-", ...
    LineWidth=1.1,MarkerSize=3,DisplayName="MDA-designed points")
yline(gramTolerance,"--",DisplayName="tolerance")
xline(mdaAssessment.gramPolicy.maximumAcceptedModeCount,":", ...
    DisplayName=sprintf("APV points accept %d",mdaAssessment.gramPolicy.maximumAcceptedModeCount))
xline(mdaDesignedAssessment.gramPolicy.maximumAcceptedModeCount,":", ...
    DisplayName=sprintf("MDA points accept %d",mdaDesignedAssessment.gramPolicy.maximumAcceptedModeCount))
if ~isempty(firstMDAGramRejection)
    plot(mdaDiagnostics.modeCount(firstMDAGramRejection), ...
        mdaDiagnostics.gramError(firstMDAGramRejection),"rx", ...
        MarkerSize=9,LineWidth=1.8,HandleVisibility="off")
end
if ~isempty(firstMDADesignedGramRejection)
    plot(mdaDesignedDiagnostics.modeCount(firstMDADesignedGramRejection), ...
        mdaDesignedDiagnostics.gramError(firstMDADesignedGramRejection),"rx", ...
        MarkerSize=9,LineWidth=1.8,HandleVisibility="off")
end
hold off
grid on
xlabel("MDA prefix mode count")
ylabel("normalized Gram error")
title("MDA G Gram error by point rule")
legend(Location="northwest")

nexttile
semilogy(apvDiagnostics.modeCount,max(apvDiagnostics.quadraticAliasingError,eps),"o-", ...
    LineWidth=1.1,MarkerSize=3,DisplayName="quadratic-product error")
hold on
yline(quadraticAliasingTolerance,"--",DisplayName="tolerance")
xline(apvAssessment.quadraticAliasingPolicy.maximumAcceptedModeCount,":", ...
    DisplayName=sprintf("accept through %d",apvAssessment.quadraticAliasingPolicy.maximumAcceptedModeCount))
if ~isempty(firstAPVQuadraticRejection)
    plot(apvDiagnostics.modeCount(firstAPVQuadraticRejection), ...
        apvDiagnostics.quadraticAliasingError(firstAPVQuadraticRejection),"rx", ...
        MarkerSize=9,LineWidth=1.8,DisplayName="first rejection")
end
hold off
grid on
xlabel("APV prefix mode count")
ylabel("normalized projection error")
title("APV coupled quadratic products")
legend(Location="northwest")

nexttile
plot(apvTransform.weights,z,"o-",LineWidth=1.1,MarkerSize=3, ...
    DisplayName="APV weights, APV points")
hold on
plot(mdaTransform.weights,z,"o-",LineWidth=1.1,MarkerSize=3, ...
    DisplayName="MDA weights, APV points")
plot(mdaDesignedTransform.weights,mdaDesignedZ,"x--",LineWidth=1.1,MarkerSize=4, ...
    DisplayName="MDA weights, MDA points")
hold off
grid on
xlabel("fitted quadrature weight (m)")
ylabel("z (m)")
title("Family metrics and point rules")
legend(Location="best")

%% Extract the retained Wave-Vortex basis arrays
% APV and MDA have different coefficient axes and different retained
% counts. APV has directly projectable F and G channels. MDA projects
% through G and uses the aligned F modes only to synthesize the diagnostic
% surface-referenced pressure.
apvFinv = apvTransform.inverseMatrix(variable="F");
apvF = apvTransform.forwardMatrix(variable="F");
apvGinv = apvTransform.inverseMatrix(variable="G");
apvG = apvTransform.forwardMatrix(variable="G");
apvModeNumber = apvTransform.modeNumber(:);
apvEquivalentDepth = apvTransform.h(:);

mdaFinv = mdaTransform.inverseMatrix(variable="F");
mdaGinv = mdaTransform.inverseMatrix(variable="G");
mdaG = mdaTransform.forwardMatrix(variable="G");
mdaModeNumber = mdaTransform.modeNumber(:);
mdaEquivalentDepth = mdaTransform.h(:);

fprintf("APV weights: sum %.6f m; minimum %.6f m.\n", ...
    sum(apvTransform.weights),min(apvTransform.weights));
fprintf("MDA weights: sum %.6f m; minimum %.6f m.\n", ...
    sum(mdaTransform.weights),min(mdaTransform.weights));
fprintf("MDA-designed weights: sum %.6f m; minimum %.6f m.\n", ...
    sum(mdaDesignedTransform.weights),min(mdaDesignedTransform.weights));
fprintf("APV retains %d modes: %d negative, zero present %d, %d positive; smallest h is %.3e m.\n", ...
    length(apvModeNumber),nnz(apvModeNumber < 0),any(apvModeNumber == 0), ...
    nnz(apvModeNumber > 0),min(apvEquivalentDepth));
fprintf("MDA retains %d modes: %d negative, zero present %d, %d positive; smallest h is %.3e m.\n", ...
    length(mdaModeNumber),nnz(mdaModeNumber < 0),any(mdaModeNumber == 0), ...
    nnz(mdaModeNumber > 0),min(mdaEquivalentDepth));
mdaDesignedModeNumber = mdaDesignedTransform.modeNumber(:);
fprintf("MDA on its own point rule retains %d modes: %d negative, zero present %d, %d positive.\n", ...
    length(mdaDesignedModeNumber),nnz(mdaDesignedModeNumber < 0), ...
    any(mdaDesignedModeNumber == 0),nnz(mdaDesignedModeNumber > 0));

%% Verify independent modal round trips
% Each retained family should reconstruct and recover its own coefficients
% to numerical precision. The MDA coefficients also synthesize an aligned
% pressure field whose value is zero at the surface by construction.
apvCoefficient = exp(-(0:length(apvModeNumber)-1).'/8) ...
    .*cos((0:length(apvModeNumber)-1).'*pi/5);
apvFRecovered = apvF*(apvFinv*apvCoefficient);
apvGRecovered = apvG*(apvGinv*apvCoefficient);
relativeAPVFError = norm(apvFRecovered-apvCoefficient)/norm(apvCoefficient);
relativeAPVGError = norm(apvGRecovered-apvCoefficient)/norm(apvCoefficient);

mdaCoefficient = exp(-(0:length(mdaModeNumber)-1).'/4) ...
    .*cos((0:length(mdaModeNumber)-1).'*pi/7);
mdaDisplacement = mdaGinv*mdaCoefficient;
mdaRecovered = mdaG*mdaDisplacement;
mdaMeanPressure = mdaFinv*mdaRecovered;
relativeMDAGError = norm(mdaRecovered-mdaCoefficient)/norm(mdaCoefficient);

fprintf("APV F coefficient round-trip error: %.3e\n",relativeAPVFError);
fprintf("APV G coefficient round-trip error: %.3e\n",relativeAPVGError);
fprintf("MDA G coefficient round-trip error: %.3e\n",relativeMDAGError);
fprintf("Maximum MDA pressure value at the surface: %.3e\n",max(abs(mdaMeanPressure(end,:))));
