%% Project partial-depth observations onto resolvable vertical modes
% Observations are an observation-specific inverse problem. The sampled
% matrix is not the full modal inverse transform; it is the modal inverse
% transform after applying the observation operator. This example uses
% irregular point samples in the upper ocean, selects a well-conditioned
% non-contiguous set of resolvable modes, and plots the true and recovered
% modal coefficients.
%
% For the G modes, the normalization
%
%     (1/g) int N2(z) eta(z)^2 dz = sum_j eta_j^2
%
% makes eta_j^2 the modal contribution to the normalized potential-energy
% spectrum. The observation operator changes what part of that spectrum is
% recoverable.

%% Build a candidate mode basis
Lz = 4000;
N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);

latitude = 31;
g = 9.81;
k = 0;
nCandidateModes = 48;
nEVP = max(256,ceil(2.1*(nCandidateModes + 1)));
zDomain = [-Lz 0];
zPlot = linspace(zDomain(1),zDomain(2),512).';

imPlot = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zPlot,latitude=latitude,nEVP=nEVP,nModes=nCandidateModes,g=g);
imPlot.normalization = Normalization.kConstant;
imPlot.upperBoundary = UpperBoundary.rigidLid;
[~,GPlot] = imPlot.modesAtWavenumber(k);

%% Define an irregular partial-depth observing grid
zObsUniform = linspace(-1600,-80,26).';
zObs = zObsUniform + 35*sin((1:length(zObsUniform)).');
zObs = sort(zObs);

imObs = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=zObs,latitude=latitude,nEVP=nEVP,nModes=nCandidateModes,g=g);
imObs.normalization = Normalization.kConstant;
imObs.upperBoundary = UpperBoundary.rigidLid;
[~,GObs] = imObs.modesAtWavenumber(k);

modeScale = max(abs(GPlot(:,1:nCandidateModes)),[],1);
PhiPlot = GPlot(:,1:nCandidateModes)./modeScale;
B = GObs(:,1:nCandidateModes)./modeScale;

%% Select modes that are identifiable by the observation operator
% The observation matrix B = H Phi maps modal coefficients to samples. QR
% with column pivoting identifies a subset of columns that is independent
% on this grid. The condition-number limit keeps the projection from
% becoming a noise amplifier.
weights = ones(size(zObs));
rankTolerance = 5e-2;
maxConditionNumber = 10;
[retainedModes,selectionDiagnostics] = selectResolvableModes(B, weights, rankTolerance, maxConditionNumber);
rejectedModes = setdiff(1:nCandidateModes, retainedModes);

[forwardProjection,resolutionMatrix,gramMatrix] = weightedProjection(B, weights, retainedModes);
retainedError = norm(resolutionMatrix(:,retainedModes) - eye(length(retainedModes)),'fro')/sqrt(length(retainedModes));
aliasingMatrix = resolutionMatrix(:,rejectedModes);
aliasingColumnNorm = vecnorm(aliasingMatrix,2,1);
maxAliasing = max(aliasingColumnNorm);

fprintf('Observation grid has %d samples over %.0f m to %.0f m.\n', length(zObs), min(zObs), max(zObs));
fprintf('Retained %d of %d candidate modes.\n', length(retainedModes), nCandidateModes);
fprintf('Retained modes: %s\n', mat2str(retainedModes));
fprintf('Condition number of retained Gram matrix: %.3e\n', cond(gramMatrix));
fprintf('Condition number after column normalization: %.3e\n', selectionDiagnostics.conditionNumber);
fprintf('Retained-mode round-trip error: %.3e\n', retainedError);
fprintf('Maximum rejected-mode aliasing norm: %.3e\n', maxAliasing);
fprintf('QR relative pivot at final retained mode: %.3e\n', selectionDiagnostics.finalRelativePivot);

%% Recover coefficients from synthetic observations
% The least-squares solve uses scaled columns for conditioning. The
% physical G-mode coefficients are the scaled coefficients divided by
% modeScale. Those physical coefficients define the potential-energy
% spectrum through eta_j^2.
trueCoefficients = zeros(nCandidateModes,1);
retainedSignalModes = retainedModes(1:min(9,length(retainedModes)));
retainedAmplitudes = [1.0 -0.65 0.42 -0.30 0.22 -0.16 0.11 -0.08 0.05].';
trueCoefficients(retainedSignalModes) = retainedAmplitudes(1:length(retainedSignalModes));

rejectedSignalModes = rejectedModes(rejectedModes <= 24);
rejectedSignalModes = rejectedSignalModes(1:min(3,length(rejectedSignalModes)));
rejectedAmplitudes = [0.12 -0.09 0.06].';
trueCoefficients(rejectedSignalModes) = rejectedAmplitudes(1:length(rejectedSignalModes));

trueProfile = PhiPlot*trueCoefficients;
observations = B*trueCoefficients;

recoveredRetainedCoefficients = forwardProjection*observations;
recoveredCoefficients = zeros(nCandidateModes,1);
recoveredCoefficients(retainedModes) = recoveredRetainedCoefficients;

predictedRecoveredCoefficients = resolutionMatrix*trueCoefficients;
unresolvedBias = resolutionMatrix(:,rejectedModes)*trueCoefficients(rejectedModes);
selectedRelativeError = norm(recoveredRetainedCoefficients - trueCoefficients(retainedModes))/norm(trueCoefficients(retainedModes));

modeScaleColumn = modeScale(:);
truePhysicalCoefficients = trueCoefficients ./ modeScaleColumn;
recoveredPhysicalCoefficients = recoveredCoefficients ./ modeScaleColumn;
truePotentialEnergySpectrum = truePhysicalCoefficients.^2;
recoveredPotentialEnergySpectrum = recoveredPhysicalCoefficients.^2;

retainedModeScale = modeScaleColumn(retainedModes);
physicalResolutionMatrix = (modeScaleColumn.' ./ retainedModeScale) .* resolutionMatrix;
spectralWindow = physicalResolutionMatrix.^2;
expectedRecoveredSpectrum = zeros(nCandidateModes,1);
expectedRecoveredSpectrum(retainedModes) = spectralWindow*truePotentialEnergySpectrum;
% The spectral window maps a true potential-energy spectrum to the expected
% recovered spectrum when modal coefficients are uncorrelated.

fprintf('\nSynthetic recovery\n');
fprintf('Relative error in retained coefficients: %.3e\n', selectedRelativeError);
fprintf('Relative unresolved-mode bias: %.3e\n', norm(unresolvedBias)/norm(trueCoefficients(retainedModes)));
fprintf('Deterministic resolution check: %.3e\n', norm(predictedRecoveredCoefficients - forwardProjection*observations));
fprintf('True normalized potential energy: %.3e\n', sum(truePotentialEnergySpectrum));
fprintf('Recovered normalized potential energy: %.3e\n', sum(recoveredPotentialEnergySpectrum));
fprintf('Expected recovered potential energy: %.3e\n', sum(expectedRecoveredSpectrum));

%% Plot true and recovered coefficients with resolution diagnostics
figure('Name','Partial-depth observational mode projection')
tiledlayout(3,2,'TileSpacing','compact','Padding','compact')

nexttile
plot(trueProfile,zPlot,'k','LineWidth',1.5), hold on
plot(observations,zObs,'ko','MarkerFaceColor','w')
xlabel('profile value')
ylabel('z (m)')
title('Synthetic partial-depth observations')
legend('truth','samples','Location','southwest')
grid on

nexttile
stem(1:nCandidateModes,truePhysicalCoefficients,'k','DisplayName','true'), hold on
stem(retainedModes,recoveredPhysicalCoefficients(retainedModes),'r','filled','DisplayName','recovered')
xlabel('mode number')
ylabel('coefficient')
title('Physical G-mode coefficients')
legend('Location','northeast')
grid on

nexttile
stem(1:nCandidateModes,truePotentialEnergySpectrum,'k','DisplayName','true'), hold on
stem(retainedModes,expectedRecoveredSpectrum(retainedModes),'b','filled','DisplayName','expected recovered')
stem(retainedModes,recoveredPotentialEnergySpectrum(retainedModes),'r','DisplayName','deterministic recovered')
xlabel('mode number')
ylabel('potential energy')
title('Resolvable potential-energy spectrum')
legend('Location','northeast')
grid on

nexttile
imagesc(1:nCandidateModes,1:length(retainedModes),abs(resolutionMatrix))
set(gca,'YDir','normal')
xlabel('true mode number')
ylabel('recovered coefficient row')
title('|A B| resolution/aliasing matrix')
colorbar

nexttile
stem(rejectedModes,aliasingColumnNorm,'filled')
xlabel('rejected mode number')
ylabel('aliasing norm')
title('Rejected-mode leakage into retained coefficients')
grid on

nexttile
imagesc(1:nCandidateModes,retainedModes,spectralWindow)
set(gca,'YDir','normal')
xlabel('true mode number')
ylabel('recovered mode number')
title('Potential-energy spectral window')
colorbar

%% Local helpers
function [retainedModes,diagnostics] = selectResolvableModes(B, weights, rankTolerance, maxConditionNumber)
weightedB = sqrt(weights) .* B;
columnNorms = vecnorm(weightedB,2,1);
normalizedB = weightedB ./ columnNorms;
[~,R,pivotOrder] = qr(normalizedB,0);
relativePivots = abs(diag(R))/abs(R(1,1));
nRank = find(relativePivots >= rankTolerance,1,'last');

selectedInPivotOrder = [];
for iPivot = 1:nRank
    candidateModes = sort([selectedInPivotOrder pivotOrder(iPivot)]);
    candidateConditionNumber = cond(normalizedB(:,candidateModes).'*normalizedB(:,candidateModes));
    if candidateConditionNumber <= maxConditionNumber
        selectedInPivotOrder = [selectedInPivotOrder pivotOrder(iPivot)];
    end
end

retainedModes = sort(selectedInPivotOrder);
diagnostics.relativePivots = relativePivots;
diagnostics.finalRelativePivot = relativePivots(length(selectedInPivotOrder));
diagnostics.conditionNumber = cond(normalizedB(:,retainedModes).'*normalizedB(:,retainedModes));
end

function [A,resolutionMatrix,gramMatrix] = weightedProjection(B, weights, retainedModes)
BS = B(:,retainedModes);
weightedBS = weights .* BS;
gramMatrix = BS.' * weightedBS;
A = gramMatrix \ weightedBS.';
resolutionMatrix = A*B;
end
