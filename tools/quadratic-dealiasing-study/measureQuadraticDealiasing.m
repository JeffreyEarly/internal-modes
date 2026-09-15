function rows = measureQuadraticDealiasing(outputPath)
% Time provider scoring alone; basis evaluation and construction are excluded.
arguments (Input)
    outputPath (1,1) string
end
maxNumCompThreads(1);
sizes = [129 32;313 64;625 128;1249 256];
rows = struct([]);
for index = 1:size(sizes,1)
    q = sizes(index,1); m = sizes(index,2);
    theta = linspace(pi,0,q).';
    F = cos(theta*(0:m-1)); G = cos(theta*(1:m));
    operation = @() assessQuadraticDealiasing(F,G,floor((q-1)/4),quadraticDealiasing="effectiveBandwidth");
    operation();
    elapsed = zeros(5,1);
    for trial = 1:numel(elapsed)
        started = tic;
        for repetition = 1:20, operation(); end
        elapsed(trial) = toc(started)/20;
    end
    row = struct(sampleCount=q,modeCount=m,batchSize=128,medianSeconds=median(elapsed),trials=elapsed,callerInputBytes=2*q*m*8);
    rows = [rows;row]; %#ok<AGROW>
    fprintf('Q=%d M=%d: %.6g seconds per call\n',q,m,row.medianSeconds);
end
writelines(jsonencode(struct(matlabVersion=string(version),threads=1,description="Two dense synthetic polynomial channels, one warmup plus five trials of twenty calls; scoring only, excluding caller basis evaluation",rows=rows),PrettyPrint=true),outputPath);
end
