function report = benchmarkBulkBasisEvaluation()
% Compare scalar, bulk, and streaming evaluation of continuous wave bases.
%
% Timings exclude construction and warmup. Array component sizes are not
% process peak memory measurements. The returned report retains raw timings
% and a numerical comparison against independent scalar evaluation.
%
% - Topic: Developer topics
% - Returns report: evaluation timings, numerical agreement, and array sizes
nEVP = 96;
nModes = 16;
nZ = 2048;
nDistinct = 24;
k = linspace(0.5e-4,5e-4,nDistinct);
requested = [k fliplr(k)];
solver = IMSolverSpectral(nEVP=nEVP);
[collection,construction] = solver.solveWaveModesAtWavenumbers(requested,N2=@(z) 1e-4*(1+0.2*exp(z/300)),zDomain=[-1000 0],f0=1e-4,nModes=nModes);
z = linspace(-1000,0,nZ).';
sampleChunkSize = 256;
pageChunkSize = 8;
peakChunkBytes = 0;
streamChecksum = 0;
streamCallbacks = 0;
% Warm every path; elapsed medians below exclude generation and warmup.
reference = scalarLoop(collection,z);
result = collection.evaluate(z,variable="F",sampleChunkSize=sampleChunkSize,pageChunkSize=pageChunkSize);
collection.evaluateChunks(z,@consume,variable="F",sampleChunkSize=sampleChunkSize,pageChunkSize=pageChunkSize);
relativeError = norm(result(:)-reference(:))/norm(reference(:));
streamRelativeChecksumError = abs(streamChecksum-sum(reference,'all'))/max(1,abs(sum(reference,'all')));
outputInfo = whos('result');
scalarTimes = zeros(1,3); evaluateTimes = scalarTimes; streamTimes = scalarTimes;
for trial = 1:3
    timer=tic; reference=scalarLoop(collection,z); scalarTimes(trial)=toc(timer);
    timer=tic; result=collection.evaluate(z,variable="F",sampleChunkSize=sampleChunkSize,pageChunkSize=pageChunkSize); evaluateTimes(trial)=toc(timer);
    streamChecksum=0; streamCallbacks=0;
    timer=tic; collection.evaluateChunks(z,@consume,variable="F",sampleChunkSize=sampleChunkSize,pageChunkSize=pageChunkSize); streamTimes(trial)=toc(timer);
end
% Selected ordinary MATLAB array allocations from the implemented kernels.
% These are component sizes, NOT an allocator/RSS peak measurement.
arrayComponents = struct('fullResultBytes',outputInfo.bytes,'streamMaximumCallbackPayloadBytes',peakChunkBytes,'streamPersistentResultBytes',0, ...
    'groupNativeColumnBytes',8*nEVP*nModes*pageChunkSize,'groupSampledValueBytes',8*sampleChunkSize*nModes*pageChunkSize, ...
    'chunkOutputBytes',8*sampleChunkSize*nModes*pageChunkSize,'oneChebyshevMatrixBytes',8*sampleChunkSize*nEVP, ...
    'physicalDerivativeDiagonalBytes',8*sampleChunkSize^2,'scalarOnePageOutputBytes',8*nZ*nModes, ...
    'scalarOneChebyshevMatrixBytes',8*nZ*nEVP,'scalarPhysicalDerivativeDiagonalBytes',8*nZ^2);
report = struct('nEVP',nEVP,'nModes',nModes,'nZ',nZ,'distinctKappaCount',nDistinct,'requestedPageCount',numel(requested), ...
    'sampleChunkSize',sampleChunkSize,'pageChunkSize',pageChunkSize,'scalarSeconds',scalarTimes,'evaluateSeconds',evaluateTimes,'streamSeconds',streamTimes, ...
    'medianScalarSeconds',median(scalarTimes),'medianEvaluateSeconds',median(evaluateTimes),'medianStreamSeconds',median(streamTimes), ...
    'relativeValueError',relativeError,'streamRelativeChecksumError',streamRelativeChecksumError,'streamCallbacks',streamCallbacks, ...
    'arrayComponents',arrayComponents,'construction',construction,'matlabVersion',version);
    function consume(values,~,~,~)
        item = whos('values');
        peakChunkBytes = max(peakChunkBytes,item.bytes);
        streamChecksum = streamChecksum+sum(values,'all');
        streamCallbacks = streamCallbacks+1;
    end
end
function result = scalarLoop(collection,z)
nColumns = numel(collection.metadata(1).columnLabels);
result = zeros(numel(z),nColumns,numel(collection.basisIndex));
for i=1:numel(collection.basisIndex)
    result(:,:,i) = collection.bases{collection.basisIndex(i)}.F(z);
end
end
