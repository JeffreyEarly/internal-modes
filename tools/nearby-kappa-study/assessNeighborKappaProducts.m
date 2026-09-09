function summary = assessNeighborKappaProducts(studyDirectory,outputDirectory,options)
% Assess unchanged neighboring bases against the same requested-kappa products.
%
% This offline helper requires the tested #19 product-inventory API on the
% MATLAB path. GG->G uses the requested wave G scalar pairing, including
% required endpoint terms. FG->F explicitly uses a plain L2 mathematical
% control, not an unavailable physical wave F projection. Original basis
% normalizations and eigenvalue-rank coordinates are retained unchanged.
% Diagnostic cross-Gram maps are reported but never applied to reuse.
%
% - Topic: Investigate nearby wavenumbers
% - Parameter studyDirectory: directory containing retained environment MAT files
% - Parameter outputDirectory: external directory for detailed assessment MAT files
% - Parameter options.caseIds: bounded retained case selection
% - Parameter options.bands: independently assessed fixed retained bands
% - Parameter options.productBudget: per case/band/channel reservation including all three outputs
% - Returns summary: compact numerical controls and explicit qualification limits
arguments
    studyDirectory (1,1) string
    outputDirectory (1,1) string
    options.caseIds (1,:) double {mustBeInteger,mustBePositive} = [1 12 17 24 29 36 52 96]
    options.bands (1,:) double {mustBeInteger,mustBePositive} = [8 16]
    options.productBudget (1,1) double {mustBeInteger,mustBePositive} = 2000
    options.sampleCount (1,1) double {mustBeInteger,mustBePositive} = 65
    options.referenceCount (1,1) double {mustBeInteger,mustBePositive} = 513
    options.quadratureCount (1,1) double {mustBeInteger,mustBePositive} = 1025
    options.referenceTolerance (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 1e-4
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
outputRecords = {};
timer = tic;
for caseId = options.caseIds
    environment = ceil(caseId/12);
    data = load(fullfile(studyDirectory,sprintf('environment-%02d.mat',environment)),'retained');
    entries = data.retained.candidateRecords;
    index = find(cellfun(@(e) e.caseId == caseId,entries),1);
    if isempty(index), error('IM21:MissingCase','Case %d is absent from the retained environment.',caseId); end
    entry = entries{index};
    requested = entry.exactBasis; neighbor = entry.sourceBasis; independent = entry.referenceBasis;
    if requested.evp.parameters.k ~= entry.record.requestedKappa || independent.evp.parameters.k ~= entry.record.requestedKappa || neighbor.evp.parameters.k ~= entry.record.sourceKappa
        error('IM21:KappaIdentity','Retained basis kappa does not match the recorded scientific problem.');
    end
    domain = requested.zDomain;
    grids = struct('id',{'sample','primary','quadrature','independent'},'z',{linspace(domain(1),domain(2),options.sampleCount).',linspace(domain(1),domain(2),options.referenceCount).',linspace(domain(1),domain(2),options.quadratureCount).',linspace(domain(1),domain(2),options.quadratureCount).'});
    for count = options.bands
        if any([numel(requested.modeNumber),numel(neighbor.modeNumber),numel(independent.modeNumber)] < count)
            error('IM21:MissingBand','Every original basis must supply the complete requested band.');
        end
        for channel = ["GGtoG","FGtoF"]
            preparationTimer = tic;
            outputVariable = "G"; inputVariables = ["G","G"];
            pairingDescription = "basis-owned wave G pairing, including required endpoints";
            if channel == "FGtoF"
                outputVariable = "F"; inputVariables = ["F","G"];
                pairingDescription = "declared plain L2 mathematical control; not a physical wave F projection";
            end
            % All deferred factors evaluate the same requested scientific
            % fields, even on the grid used by the independent output solve.
            factors = cell(1,2);
            for side = 1:2
                variable = inputVariables(side);
                factors{side} = struct('id',"source"+side,'family',"requestedWave"+variable,'labels',string(requested.modeNumber(1:count)),'ordinals',1:count,'frequencySigns',zeros(1,count),'countRole',"retained",'evaluate',@(z,columns,gridId) sourceValues(requested,variable,z,columns,gridId),'provenance',struct('requestedKappa',entry.record.requestedKappa,'normalization',string(requested.normalization),'variable',variable,'sourceSubstitution',false));
            end
            labels = string(requested.modeNumber(1:count));
            status = "qualified";
            if ~entry.record.referenceQualifiedExperimental, status = "inconclusive"; end
            % Each plan assesses one fixed band. The difference dual is the
            % difference of these complete-band projection matrices, not a
            % prefix of a larger band's difference matrix.
            outputs = cell(1,3);
            for method = 1:3
                outputs{method} = struct('id',string(method),'family',"requestedCoordinates"+outputVariable,'labels',labels,'ordinals',1:count,'countRole',"retained",'prepare',@(g) prepareOutput(g,requested,neighbor,independent,count,outputVariable,method,status),'provenance',struct('method',method,'pairing',pairingDescription,'coordinateMapping',"eigenvalue rank; unchanged original normalization/orientation"));
            end
            products = table(repmat("case"+caseId,3,1),repmat(channel,3,1),ones(3,1),2*ones(3,1),(1:3).',VariableNames=["interactionId","channel","factorA","factorB","output"]);
            inventory = IMProductInventory(factors,products,outputs);
            plan = inventory.fixedPlan(prefixCounts=count,productBudget=options.productBudget);
            metadataPreparationSeconds = toc(preparationTimer);
            assessment = plan.assess(grids,chunkSize=64,minimumReciprocalCondition=1e-12);
            e = assessment.evidence;
            finite = isfinite(e{1}.error) & isfinite(e{2}.error) & isfinite(e{3}.error);
            slack = 1e-10*max(1,max([e{1}.error(finite),e{2}.error(finite),e{3}.error(finite)],[],"all"));
            if any(e{2}.error(finite) > e{1}.error(finite)+e{3}.error(finite)+slack)
                error('IM21:CoefficientTriangle','Baseline, reuse, and total coefficient errors do not share one norm/reference.');
            end
            if entry.record.sourceKappa == entry.record.requestedKappa && any(e{3}.error(finite) ~= 0)
                error('IM21:ExactRepeatReuse','An exact-repeat basis must have zero reuse coefficient difference.');
            end
            baseline = max(e{1}.error,[],'all'); total = max(e{2}.error,[],'all'); reuse = max(e{3}.error,[],'all');
            referenceError = max([e{1}.referenceError,e{2}.referenceError],[],'all');
            quadratureError = max([e{1}.quadratureError,e{2}.quadratureError],[],'all');
            independentError = max([e{1}.independentSolveError,e{2}.independentSolveError],[],'all');
            referenceQualified = status == "qualified" && referenceError <= options.referenceTolerance;
            diagnosticsTimer = tic;
            transform = transformDiagnostics(requested,neighbor,count,outputVariable,grids(3).z);
            transformDiagnosticsSeconds = toc(diagnosticsTimer);
            row = struct('caseId',caseId,'profile',string(entry.record.profile),'surface',string(entry.record.surface),'sourceKappa',entry.record.sourceKappa,'requestedKappa',entry.record.requestedKappa,'relativeOffset',entry.record.relativeOffset,'count',count,'channel',channel,'pairing',pairingDescription,'baselineAliasing',baseline,'reuseCoefficientDifference',reuse,'neighborTotalError',total,'referenceError',referenceError,'quadratureError',quadratureError,'independentSolveError',independentError,'referenceQualifiedExperimental',referenceQualified,'referenceTolerance',options.referenceTolerance,'reservedProducts',plan.reservedProducts,'evaluatedProducts',assessment.costs.evaluatedProducts,'examinedProductsPerOutput',numel(e{1}.isZero),'assessmentSeconds',assessment.costs.totalSeconds,'metadataPreparationSeconds',metadataPreparationSeconds,'transformDiagnosticsSeconds',transformDiagnosticsSeconds,'synthesisDifference',transform.synthesisDifference,'outsideRequestedSpan',transform.outsideRequestedSpan,'unchangedOrientationMatched',all(transform.overlapSigns==1),'requestedLabels',string(requested.modeNumber(1:count)),'neighborLabels',string(neighbor.modeNumber(1:count)),'independentLabels',string(independent.modeNumber(1:count)),'requestedNormalization',string(requested.normalization),'neighborNormalization',string(neighbor.normalization),'independentNormalization',string(independent.normalization),'coordinateMapping',"eigenvalue rank; no rotation, fit, rescaling, or sign change applied",'productionQualification',false);
            if ~row.unchangedOrientationMatched || row.requestedNormalization ~= row.neighborNormalization || row.requestedNormalization ~= row.independentNormalization
                row.referenceQualifiedExperimental = false;
            end
            outputRecords{end+1} = row; %#ok<AGROW>
            save(fullfile(outputDirectory,sprintf('case-%03d-band-%02d-%s.mat',caseId,count,channel)),'assessment','row','transform','-v7.3');
            fprintf('products case=%d band=%d %s base=%g reuse=%g total=%g reference=%g qualified=%d\n',caseId,count,channel,baseline,reuse,total,referenceError,row.referenceQualifiedExperimental);
        end
    end
end
rows = [outputRecords{:}];
summary = struct('caseIds',options.caseIds,'bands',options.bands,'completedControls',numel(rows),'qualifiedControls',nnz([rows.referenceQualifiedExperimental]),'elapsedSeconds',toc(timer),'options',options,'scope',"Bounded generic GG-to-G and plain-L2 FG-to-F numerical controls; no WVM physical source inventory, Fourier interaction coverage, superposition bound, or pressure/momentum qualification",'productionQualification',false);
save(fullfile(outputDirectory,'product-summary.mat'),'summary','rows');
fid=fopen(fullfile(outputDirectory,'product-summary.json'),'w'); cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(struct('summary',summary,'rows',rows),PrettyPrint=true));
end

function values = sourceValues(basis,variable,z,columns,~)
values = basis.(variable)(z);
values = values(:,columns);
end

function context = prepareOutput(grids,requested,neighbor,independent,count,variable,method,status)
labels = string(requested.modeNumber(1:count));
[Ar,Wr] = observed(requested,variable,grids(2).z,count);
primaryTarget = symmetric(Ar.'*Wr*Ar);
majorant = primaryTarget; % Both declared controls have positive metrics.
[A,W] = observed(requested,variable,grids(1).z,count);
[B,Wb] = observed(neighbor,variable,grids(1).z,count);
[Bref,Wbr] = observed(neighbor,variable,grids(2).z,count);
exact = IMProjection(A,full(W),primaryTarget,majorantGramMatrix=majorant,columnLabels=labels);
reused = IMProjection(B,full(Wb),symmetric(Bref.'*Wbr*Bref),majorantGramMatrix=majorant,columnLabels=labels);
if method == 1
    projection = exact;
elseif method == 2
    projection = reused;
else
    projection = IMProjection.fromPairing(reused.forwardMatrix-exact.forwardMatrix,eye(count),eye(count),majorantGramMatrix=majorant,columnLabels=labels,provenance=struct('kind',"unchanged-neighbor minus requested fixed-band coefficient operator",'coordinateMapApplied',false));
end
references = cell(1,3);
roles = ["primary","quadrature","independent"];
for i = 1:3
    z = grids(i+1).z;
    outputBasis = requested;
    if i == 3, outputBasis = independent; end
    [C,metric] = observed(outputBasis,variable,z,count);
    % Product norms always use the requested physical observation metric,
    % independently of which solved output basis supplies the signed dual.
    [~,normMetric] = observed(requested,variable,z,count);
    pairing = C.'*metric;
    target = symmetric(C.'*metric*C);
    referenceMajorant = target;
    if method == 3
        pairing = zeros(count,numel(z)); target = eye(count); referenceMajorant = majorant;
    end
    references{i} = struct('gridId',string(grids(i+1).id),'pairingMatrix',pairing,'normMatrix',normMetric,'targetGramMatrix',target,'majorantGramMatrix',referenceMajorant,'role',roles(i),'status',status,'provenance',struct('outputKappa',outputBasis.evp.parameters.k,'originalLabels',outputBasis.modeNumber(1:count),'normalization',string(outputBasis.normalization),'independentSolve',i==3,'coordinateMapping',"eigenvalue rank, original normalized/oriented columns",'referenceGridCount',numel(z),'differenceDualReferenceIsZero',method==3));
end
context = struct('projection',projection,'references',{references});
end

function [values,metric] = observed(basis,variable,z,count)
values = basis.(variable)(z); values = values(:,1:count);
weights = [(z(2)-z(1))/2;(z(3:end)-z(1:end-2))/2;(z(end)-z(end-1))/2];
if variable == "G"
    recipe = basis.projectionRecipe(variable="G");
    metric = sparse(recipe.metric(z,weights));
else
    metric = spdiags(weights,0,numel(z),numel(z));
end
end

function diagnostics = transformDiagnostics(requested,neighbor,count,variable,z)
[A,W] = observed(requested,variable,z,count);
[B,~] = observed(neighbor,variable,z,count);
majorant = symmetric(A.'*W*A);
R = chol(majorant);
difference = (B-A)/R;
synthesisDifference = sqrt(max(0,max(eig(symmetric(difference.'*W*difference)))));
crossGramMap = majorant\(A.'*W*B);
outside = (B-A*crossGramMap)/R;
outsideRequestedSpan = sqrt(max(0,max(eig(symmetric(outside.'*W*outside)))));
diagnostics = struct('synthesisDifference',synthesisDifference,'outsideRequestedSpan',outsideRequestedSpan,'crossGramMap',crossGramMap,'crossGramMapAppliedToReuse',false,'overlapSigns',sign(diag(A.'*W*B)).','norm',"positive requested observation metric and requested coefficient metric",'sampledGramRequested',majorant,'sampledGramNeighbor',symmetric(B.'*W*B));
end

function matrix = symmetric(matrix)
matrix = (matrix+matrix.')/2;
end
