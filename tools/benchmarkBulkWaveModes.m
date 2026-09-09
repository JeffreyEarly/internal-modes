function report = benchmarkBulkWaveModes(options)
% Compare independent scalar and exact bulk continuous-wave construction.
%
% Scalar baselines time all requested pages and only distinct wavenumbers
% separately, distinguishing exact deduplication from shared preparation.
% Bulk stage timings come from the actual measured solve. Shared assembly
% storage is an allocation estimate, not process peak memory. Evaluation
% timings and process peak RSS belong to the separate collection benchmark.
% No benchmark timing is an accuracy or model-qualification claim.
%
% - Topic: Developer topics
% - Parameter options.nEVP: native spectral resolution
% - Parameter options.nModes: columns per positive-wavenumber page
% - Parameter options.repetitions: measured repetitions after warmup
% - Parameter options.coordinateKind: spectral coordinate
% - Returns report: raw construction timings and exact mapping counts
arguments (Input)
    options.nEVP (1,1) double {mustBeInteger,mustBePositive} = 64
    options.nModes (1,1) double {mustBeInteger,mustBePositive} = 8
    options.repetitions (1,1) double {mustBeInteger,mustBePositive} = 3
    options.coordinateKind (1,1) string {mustBeMember(options.coordinateKind,["z","wkb","density"])} = "z"
end
arguments (Output)
    report (1,1) struct
end
N2 = @(z) 1e-4*(1+0.2*exp(z/300));
distinctKappa = linspace(1e-4,1e-3,8);
kappa = repmat(distinctKappa,1,2);
solver = IMSolverSpectral(nEVP=options.nEVP,coordinateKind=options.coordinateKind);
factory = @(k) IMInternalModes.waveModesAtWavenumber(k=k,N2=N2,zDomain=[-1000 0],f0=1e-4);
solver.solveEVP(factory(kappa(1)),nModes=options.nModes);
solver.solveWaveModesAtWavenumbers(kappa,N2=N2,zDomain=[-1000 0],f0=1e-4,nModes=options.nModes);
scalarSeconds = zeros(options.repetitions,1);
scalarDistinctSeconds = scalarSeconds;
bulkSeconds = scalarSeconds;
bulkDiagnostics = cell(options.repetitions,1);
for i = 1:options.repetitions
    timer = tic;
    for k = kappa
        solver.solveEVP(factory(k),nModes=options.nModes);
    end
    scalarSeconds(i) = toc(timer);
    timer = tic;
    for k = distinctKappa
        solver.solveEVP(factory(k),nModes=options.nModes);
    end
    scalarDistinctSeconds(i) = toc(timer);
    timer = tic;
    [~,bulkDiagnostics{i}] = solver.solveWaveModesAtWavenumbers(kappa,N2=N2,zDomain=[-1000 0],f0=1e-4,nModes=options.nModes);
    bulkSeconds(i) = toc(timer);
end
report = struct("nEVP",options.nEVP,"nModes",options.nModes,"coordinateKind",options.coordinateKind,"kappa",kappa, ...
    "scalarSeconds",scalarSeconds,"scalarDistinctSeconds",scalarDistinctSeconds,"bulkSeconds",bulkSeconds, ...
    "medianRequestedSpeedup",median(scalarSeconds)/median(bulkSeconds),"medianDistinctSpeedup",median(scalarDistinctSeconds)/median(bulkSeconds), ...
    "bulkDiagnostics",{bulkDiagnostics},"processPeakMemory","not measured","evaluationTime","not measured");
report.solveCounts = struct("scalarRequested",numel(kappa),"scalarDistinct",numel(distinctKappa),"bulk",numel(distinctKappa));
stages = ["setupSeconds","sharedAssemblySeconds","assemblySeconds","eigensolveSeconds","finalizationSeconds","collectionSeconds","totalSeconds"];
medianSeconds = zeros(numel(stages),1);
for iStage = 1:numel(stages)
    stage = stages(iStage);
    medianSeconds(iStage) = median(cellfun(@(entry) entry.(stage),bulkDiagnostics));
end
report.bulkStageMedians = table(stages(:),medianSeconds,VariableNames=["stage","seconds"]);
end
