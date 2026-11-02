function [measurements,environment] = measureSharedRecovery()
% Measure prepared recovery setup, cached application and value payload.
%
% Run with the authoring repository on the path and one MATLAB thread. The
% synthetic complex fixtures compare construction policies, not physical
% grid quality or maximum usable bands. No EVP or weight fit is timed.
%
% - Returns measurements: separate setup/application seconds and payload bytes
% - Returns environment: runtime, platform and fixture provenance
arguments (Output)
    measurements table
    environment struct
end

previousRandomState = rng;
cleanup = onCleanup(@() rng(previousRandomState));
rng(57,"twister");
sampleCounts = [65 195];
columnCounts = [24 96];
policies = ["sampledLeastSquares" "correctedPhysicalDual" "extendedSpan"];
rows = cell(6,10);
row = 0;
batchCount = 256;
for sizeIndex = 1:numel(sampleCounts)
    nSamples = sampleCounts(sizeIndex);
    nColumns = columnCounts(sizeIndex);
    nExtended = ceil(1.5*nColumns);
    [frame,~] = qr(randn(nSamples)+1i*randn(nSamples));
    amplitudes = 10.^linspace(-2,2,nExtended);
    synthesis = frame(:,1:nColumns).*amplitudes(1:nColumns);
    physicalDual = (frame(:,1:nColumns)'+0.1*frame(:,nColumns+(1:nColumns))')./amplitudes(1:nColumns).';
    labels = "mode:"+string(1:nColumns);
    coefficients = randn(nColumns,batchCount)+1i*randn(nColumns,batchCount);
    values = randn(nSamples,batchCount)+1i*randn(nSamples,batchCount);
    for policy = policies
        options = struct("policy",policy,"columnLabels",labels, ...
            "coefficientScale",1./amplitudes(1:nColumns), ...
            "observationFactor",eye(nSamples));
        if policy == "correctedPhysicalDual"
            options.physicalDual = physicalDual;
        elseif policy == "extendedSpan"
            options.extendedSynthesis = frame(:,1:nExtended).*amplitudes;
            options.extendedColumnLabels = "mode:"+string(1:nExtended);
            options.extendedActiveColumnMask = true(1,nExtended);
            options.extendedCoefficientScale = 1./amplitudes;
            options.retainedColumnIndices = 1:nColumns;
        end
        pairs = namedargs2cell(options);
        build = @() IMDiscreteTransform.fromObservationMatrices(synthesis,pairs{:});
        transform = build();
        setupSeconds = timeit(build);
        forwardSeconds = timeit(@() transform.transformForward(values));
        backSeconds = timeit(@() transform.transformBack(coefficients));
        payload = struct("synthesis",transform.synthesisMatrix, ...
            "forward",transform.forwardMatrix,"physicalDual",transform.physicalDualMatrix, ...
            "specification",transform.recoverySpecification,"diagnostics",transform.diagnostics); %#ok<NASGU> Used by whos below.
        payloadInfo = whos("payload");
        row = row+1;
        rows(row,:) = {policy,nSamples,nColumns,nExtended,batchCount, ...
            setupSeconds,forwardSeconds,backSeconds,payloadInfo.bytes, ...
            transform.diagnostics.roundTripError};
    end
end
measurements = cell2table(rows,"VariableNames",["policy","sampleCount","columnCount", ...
    "fixtureExtendedCount","batchCount","setupSeconds","forwardBatchSeconds", ...
    "backBatchSeconds","valuePayloadBytes","roundTripError"]);
environment = struct("matlabVersion",string(version),"platform",string(computer), ...
    "threadCount",maxNumCompThreads,"seed",57,"batchCount",batchCount, ...
    "memoryMeaning","whos bytes of explicit matrix/specification/diagnostic struct; excludes object/runtime overhead and is not process RSS");
end
