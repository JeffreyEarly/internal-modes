function report = runModeOrientationStudy(scientificResultPath,outputJSON)
% Time orientation alone on the pinned 1,140 already-solved wave bases.
% Configure the intended provider and released dependencies before calling.
arguments (Input)
    scientificResultPath (1,1) string
    outputJSON (1,1) string
end
loaded = load(scientificResultPath,"assessment");
kappa = [0 reshape(loaded.assessment.pages.kappa,1,[])];
assert(numel(kappa)==570,"Expected the pinned 570-page construction workload.");
N2 = @(z) 1e-4*exp(2*z/700);
surface = IMBoundaryCondition(a=0,b=1,c=1,d=0);
bases = cell(1,2);
resolutions = [104 156];
for index = 1:2
    collection = IMSolverSpectral(nEVP=resolutions(index),coordinateKind="wkb").solveWaveModesAtWavenumbers(kappa, ...
        N2=N2,zDomain=[-1000 0],f0=2*7.2921e-5*sind(30),g=9.81,surfaceBoundary=surface,nModes=64,nInertialModes=64);
    bases{index} = collection.bases;
end
bases = [bases{:}];
expected = orientAll(bases);
seconds = zeros(3,1);
checksums = zeros(3,1);
for repetition = 1:3
    started = tic;
    checksums(repetition) = orientAll(bases);
    seconds(repetition) = toc(started);
    fprintf('Orientation repetition %d: %.6f s\n',repetition,seconds(repetition));
end
assert(all(checksums==expected),"Orientation checksum changed between trials.");
report = struct(matlabVersion=version,providerRoot=fileparts(fileparts(which('IMInternalModesBasis'))), ...
    nBases=numel(bases),resolutions=resolutions,nModes=64,seconds=seconds,medianSeconds=median(seconds),checksums=checksums, ...
    scope="One warmup and three orientation-only trials over already-solved, already-oriented bases. Includes raw G/F evaluation, derivative, policy, sign application, metadata and checksum; excludes eigensolves, collection construction and full WVM construction.");
writelines(jsonencode(report,PrettyPrint=true),outputJSON);
end

function checksum = orientAll(bases)
checksum = 0;
for index = 1:numel(bases)
    oriented = bases{index}.orientModeSigns();
    checksum = checksum+sum(oriented.nativeModes,"all");
end
end
