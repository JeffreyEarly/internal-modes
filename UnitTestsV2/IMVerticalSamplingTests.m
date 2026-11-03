classdef IMVerticalSamplingTests < matlab.unittest.TestCase
    % Independent references for caller-prescribed scalar observations.

    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function pointMeansAndIntegralsHaveIndependentFourierReferences(testCase)
            harmonics = [0 1 3 7];
            basis = fourierBasis(harmonics);
            points = IMVerticalSampling.pointValues("F",[-1.81;-.22],observationLabels=["F:lower","F:upper"],units="1");
            means = IMVerticalSampling.cellAverages("F",[-1.43 -.63],observationLabels="F:mean",units="1");
            integrals = IMVerticalSampling.cellIntegrals("G",[-1.71 -.31],observationLabels="G:integral",units="m");
            derivative = IMVerticalSampling.pointValues("Fz",-.73,observationLabels="Fz:point",units="m^-1");
            observations = [means.observations;points.observations;integrals.observations;derivative.observations];
            sampling = IMVerticalSampling(observations);
            [A,record] = sampling.evaluate(basis);
            expected = exactFourierRows(observations,harmonics,basis.zDomain);
            testCase.verifyEqual(A,expected,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(record.observations,observations);
            testCase.verifyEqual(record.observationLabels,observations.observationLabel.');
            testCase.verifyEqual(record.columnLabels,basis.columnLabels);
            testCase.verifyEqual(record.resolvedVariables(end),"F");
            testCase.verifyEqual(record.resolvedDerivativeOrders(end),1);
            midpoint = basis.F(mean(observations{1,["zLower","zUpper"]}));
            testCase.verifyGreaterThan(norm(A(1,:)-midpoint),.1);
            testCase.verifyEqual(A(1,1),1,AbsTol=3e-14);
            testCase.verifyEqual(A(4,1),0,AbsTol=0);
            meanG = IMVerticalSampling.cellAverages("G",[-1.71 -.31],observationLabels="G:mean");
            testCase.verifyEqual(A(4,:),1.4*meanG.evaluate(basis),RelTol=3e-12,AbsTol=3e-12);
            renamed = observations([4 2 1 2 5],:);
            renamed.observationLabel = "request:"+string((1:height(renamed)).');
            reordered = IMVerticalSampling(renamed);
            testCase.verifyEqual(reordered.evaluate(basis),expected([4 2 1 2 5],:),RelTol=3e-12,AbsTol=3e-12);
        end

        function selectedHighHarmonicNeedsResolvedIntegrationAndFiniteBudget(testCase)
            harmonics = [0 1 2 3 137];
            full = fourierBasis(harmonics);
            basis = full.selectColumns(5);
            means = IMVerticalSampling.cellAverages("F",[-1.774 -.346],observationLabels="high:mean");
            integrals = IMVerticalSampling.cellIntegrals("G",[-1.653 -.137],observationLabels="high:integral");
            sampling = IMVerticalSampling([means.observations;integrals.observations]);
            [A,record] = sampling.evaluate(basis);
            expected = exactFourierRows(sampling.observations,137,basis.zDomain);
            testCase.verifyEqual(A,expected,AbsTol=3e-11);
            testCase.verifyEqual(record.columnLabels,full.columnLabels(5));
            testCase.verifyTrue(all(isfinite(record.integration.estimatedError),"all"));
            testCase.verifyTrue(all(record.integration.finalOrder <= record.integration.maxOrder));
            testCase.verifyError(@() sampling.evaluate(basis,initialOrder=4,maxOrder=8,relativeTolerance=1e-13),"IMVerticalSampling:IntegrationBudgetExceeded");
            testCase.verifyError(@() sampling.evaluate(basis,initialOrder=8,maxOrder=8),"IMVerticalSampling:IntegrationBudgetExceeded");
        end

        function staggeredChannelsRecoverSharedComplexCoordinates(testCase)
            basis = fourierBasis(0:3);
            observations = prescribedFixture();
            sampling = IMVerticalSampling(observations);
            [A,record] = sampling.evaluate(basis);
            testCase.verifyEqual(A,exactFourierRows(observations,0:3,basis.zDomain),RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(A(observations.variable == "G",1),zeros(3,1),AbsTol=0);
            coefficients = [1+.7i;-.4+.2i;.3-1.1i;1.2-.6i];
            recovery = IMDiscreteTransform.fromObservationMatrices(A,columnLabels=record.columnLabels,activeColumnMask=true(1,4),fitFactor=eye(size(A,1)),metadata=record);
            testCase.verifyTrue(recovery.success);
            testCase.verifyEqual(recovery.activeColumnMask,true(1,4));
            testCase.verifyEqual(recovery.transformForward(A*coefficients),coefficients,AbsTol=3e-12);
            testCase.verifyEqual(recovery.transformBack(coefficients),A*coefficients,AbsTol=0);
            testCase.verifyEmpty(recovery.physicalDualMatrix);
            testCase.verifyTrue(isnan(recovery.diagnostics.observationNoiseGain));
            extraG = IMVerticalSampling.cellAverages("G",[-1.47 -.23],observationLabels="G:fourth");
            onlyG = IMVerticalSampling([observations(observations.variable == "G",:);extraG.observations]);
            gA = onlyG.evaluate(basis);
            testCase.verifySize(gA,[4 4]);
            testCase.verifyEqual(rank(gA),3);
            unresolved = IMDiscreteTransform.fromObservationMatrices(gA,columnLabels=basis.columnLabels,activeColumnMask=true(1,4),fitFactor=eye(size(gA,1)),failureMode="diagnostic");
            testCase.verifyFalse(unresolved.success);
            testCase.verifyEqual(unresolved.activeColumnMask,true(1,4));
            testCase.verifyEmpty(unresolved.forwardMatrix);
            testCase.verifyTrue(all(observations.zLower > basis.zDomain(1) & observations.zUpper < basis.zDomain(2)));
        end

        function neutralCellsSplitSupportsAndKeepClassicalDerivativeSemantics(testCase)
            intervals = [-1.75 -1.1;-.65 -.05];
            background = IMModeBackground(N2=@(z) double(z < -1.75 | (z > -1.1 & z < -.65) | z > -.05),zDomain=[-2 0],g=2,rotationRate=.6/(2*sind(33)),latitude=33);
            basis = IMModes.neutralLabelModes(background,neutralIntervals=intervals,nModesPerInterval=[2 3],normalization="maxAbsG");
            means = IMVerticalSampling.cellAverages("G",[-1.93 -.3;-1.6 -.2;-.9 -.03],observationLabels=["G:cross:1","G:cross:2","G:cross:3"]);
            classical = IMVerticalSampling.cellIntegrals("G",[-2 0],derivativeOrder=2,observationLabels="Gzz:classical");
            apv = IMVerticalSampling.cellIntegrals("APV",[-1.6 -.2],observationLabels="APV:cross");
            sampling = IMVerticalSampling([means.observations;classical.observations;apv.observations]);
            [A,record] = sampling.evaluate(basis,initialOrder=16,maxOrder=256);
            expected = exactNeutralRows(basis,sampling.observations,intervals);
            testCase.verifyEqual(A,expected,RelTol=5e-11,AbsTol=5e-11);
            testCase.verifyGreaterThan(norm(A(4,:)),1);
            testCase.verifyEqual(record.derivativeConvention,"classicalIntervalInterior");
            testCase.verifyGreaterThanOrEqual(record.integration.subintervalCount(4),5);
            testCase.verifyEqual(A(5,:),-background.f0*exactNeutralRows(basis,replaceWithGz(apv.observations),intervals),RelTol=5e-11,AbsTol=5e-11);
            gap = IMVerticalSampling.cellAverages("G",[-1.08 -.68],observationLabels="G:gap");
            testCase.verifyEqual(gap.evaluate(basis),zeros(1,basis.nColumns),AbsTol=0);
            edge = intervals(1,1);
            thin = IMVerticalSampling.cellAverages("APV",[edge-2*eps(edge),edge],observationLabels="APV:thin:outside");
            testCase.verifyError(@() thin.evaluate(basis),"IMVerticalSampling:InvalidSupport");
        end

        function finiteDifferenceCompanionCellsMatchExactPchipIntegrals(testCase)
            zF = -2+2*linspace(0,1,19).'.^1.5;
            zG = [zF(1);.5*(zF(1:end-1)+zF(2:end));zF(end)];
            background = IMModeBackground.constantStratification(N0=1,zDomain=[-2 0],g=2,rotationRate=.6/(2*sind(33)),latitude=33);
            basis = IMModes.geostrophicAPVModes(background,energyBoundary=IMGeneralizedEnergyBoundary(g0=-.4,gd=.3),nPositiveModes=2,normalization="maxAbsF",solver=IMSolverFiniteDifference(z=zF));
            f = IMVerticalSampling.cellIntegrals("F",[-1.803 -.173],observationLabels="F:FD");
            g = IMVerticalSampling.cellIntegrals("G",[-1.873 -.237;-1.17 -.031],observationLabels=["G:FD:1","G:FD:2"]);
            sampling = IMVerticalSampling([g.observations;f.observations]);
            [A,record] = sampling.evaluate(basis,initialOrder=4,maxOrder=16);
            expected = [integratePchip(zG,basis.G(zG),[-1.873 -.237]);integratePchip(zG,basis.G(zG),[-1.17 -.031]);integratePchip(zF,basis.F(zF),[-1.803 -.173])];
            testCase.verifyEqual(A,expected,RelTol=3e-11,AbsTol=3e-12);
            requiredPieces = 1+sum(zG > -1.873 & zG < -.237);
            testCase.verifyGreaterThanOrEqual(record.integration.subintervalCount(1),requiredPieces);
        end

        function normalizationRotationSelectionAndPageRecordsKeepTheirMeaning(testCase)
            observations = prescribedFixture();
            sampling = IMVerticalSampling(observations);
            background = IMModeBackground.constantStratification(N0=1,zDomain=[-2 0],g=2,rotationRate=.6/(2*sind(33)),latitude=33);
            original = IMModes.geostrophicZeroAPVModesAtWavenumber(background,.7,energyBoundary=IMGeneralizedEnergyBoundary(g0=0,gd=0),solver=IMSolverAnalytical());
            rotated = original.rotateWithPencil(name="sampling:rotation",leftMatrix=[3 .4;.4 1],rightMatrix=[2 .2;.2 1]);
            changed = rotated.withNormalization("maxAbsG");
            T = original.sourceCoordinateMap\changed.sourceCoordinateMap;
            [A,record] = sampling.evaluate(original);
            [newA,newRecord] = sampling.evaluate(changed);
            testCase.verifyEqual(newA,A*T,RelTol=3e-11,AbsTol=3e-12);
            testCase.verifyEqual(record.columnLabels,original.columnLabels);
            testCase.verifyEqual(newRecord.columnLabels,changed.columnLabels);
            testCase.verifyEqual(newRecord.scientificSpecification.resolved.sourceCoordinateMap,changed.sourceCoordinateMap,AbsTol=0);
            selected = changed.selectColumns(2);
            [selectedA,selectedRecord] = sampling.evaluate(selected);
            testCase.verifyEqual(selectedA,newA(:,2),RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(selectedRecord.columnLabels,changed.columnLabels(2));
            first = IMModes.waveModesAtWavenumber(background,.3,nInternalModes=2,solver=IMSolverAnalytical());
            second = IMModes.waveModesAtWavenumber(background,.8,nInternalModes=2,solver=IMSolverAnalytical());
            collection = IMBasisCollection({first,second},basisIndex=[2 1 2]);
            values = cell(1,3); records = cell(1,3);
            for page = [3 1 2]
                basis = collection.bases{collection.basisIndex(page)};
                [values{page},records{page}] = sampling.evaluate(basis);
                testCase.verifyEqual(records{page}.scientificSpecification,basis.scientificSpecification);
                testCase.verifyEqual(records{page}.constructionRecord,basis.constructionRecord);
                testCase.verifyEqual(records{page}.columnLabels,basis.columnLabels);
            end
            testCase.verifyEqual(values{1},values{3},AbsTol=0);
            testCase.verifyEqual(records{1}.scientificSpecification.resolved.horizontal.k,.8);
            testCase.verifyEqual(records{2}.scientificSpecification.resolved.horizontal.k,.3);
            testCase.verifyGreaterThan(norm(values{1}-values{2},"fro"),.01);
        end

        function unseenFieldComponentsAreNotCertifiedByObservationRank(testCase)
            full = fourierBasis(0:11);
            basis = full.selectColumns(1:4);
            sampling = IMVerticalSampling(prescribedFixture());
            fullA = sampling.evaluate(full);
            [A,record] = sampling.evaluate(basis);
            recovery = IMDiscreteTransform.fromObservationMatrices(A,columnLabels=record.columnLabels,activeColumnMask=true(1,4),fitFactor=eye(size(A,1)),metadata=record);
            testCase.verifyTrue(recovery.success);
            hidden = full.F(-.739);
            invisible = null(fullA);
            direction = invisible*(invisible'*hidden');
            direction = direction/norm(direction);
            coefficients = [1+.3i;-.7+.2i;.2-.6i;.8+.9i];
            fullCoefficients = [coefficients;zeros(full.nColumns-4,1)]+(.7+.4i)*direction;
            source = fullA*fullCoefficients;
            testCase.verifyEqual(source,A*coefficients,AbsTol=3e-13);
            recovered = recovery.transformForward(source);
            testCase.verifyEqual(recovered,coefficients,AbsTol=3e-12);
            testCase.verifyGreaterThan(abs(hidden*fullCoefficients-basis.F(-.739)*recovered),.1);
        end

        function emptySchemasAndImmutableSnapshotsDoNotInventObservations(testCase)
            basis = fourierBasis(0:3);
            emptyRows = IMVerticalSampling.pointValues("F",zeros(0,1));
            [A,record] = emptyRows.evaluate(basis);
            testCase.verifySize(A,[0 4]);
            testCase.verifyEmpty(record.observationLabels);
            sampling = IMVerticalSampling(prescribedFixture());
            [emptyColumns,emptyRecord] = sampling.evaluate(basis.selectColumns([]));
            testCase.verifySize(emptyColumns,[6 0]);
            testCase.verifyEmpty(emptyRecord.columnLabels);
            bad = sampling.observations;
            bad.observationLabel(2) = bad.observationLabel(1);
            testCase.verifyError(@() IMVerticalSampling(bad),"IMVerticalSampling:InvalidObservationLabels");
            bad = sampling.observations;
            bad.zUpper(1) = bad.zLower(1)+.1;
            testCase.verifyError(@() IMVerticalSampling(bad),"IMVerticalSampling:InvalidSupport");
            testCase.verifyError(@() IMVerticalSampling.cellAverages("F",[-1 -1+eps(1)/2],observationLabels="tiny"),"IMVerticalSampling:InvalidSupport");
            outside = IMVerticalSampling.pointValues("F",-2-eps(2));
            testCase.verifyError(@() outside.evaluate(basis),"IMVerticalSampling:OutsideDomain");
            unavailable = IMVerticalSampling.pointValues("APV",-1);
            testCase.verifyError(@() unavailable.evaluate(basis),"IMVerticalSampling:UnavailableVariable");
            unsupported = IMVerticalSampling.pointValues("F",-1,derivativeOrder=2);
            testCase.verifyError(@() unsupported.evaluate(basis),"IMVerticalSampling:UnsupportedDerivative");
            counter = IMScalarTransformTestCounter();
            background = IMModeBackground.constantStratification(N0=1,zDomain=[-2 0],g=2,rotationRate=.6/(2*sind(33)),latitude=33,surfaceCondition="rigidLid");
            source = IMModes.hydrostaticFModes(background,nInternalModes=3,solver=IMScalarTransformCountingSolver(counter));
            initialCount = counter.solveCount;
            [snapshot,identity] = sampling.evaluate(source);
            recovery = IMDiscreteTransform.fromObservationMatrices(snapshot,columnLabels=identity.columnLabels,activeColumnMask=true(1,4),fitFactor=eye(size(snapshot,1)),metadata=identity);
            testCase.verifyTrue(recovery.success);
            coefficients = [1+.4i;-.3+.9i;.7-.6i;-.2+.1i];
            expected = snapshot*coefficients;
            altered = source.withNormalization("maxAbsF").selectColumns([4 2]);
            testCase.verifyEqual(altered.nColumns,2);
            for j = 1:3
                testCase.verifyEqual(recovery.transformBack(coefficients),expected,AbsTol=0);
                testCase.verifyEqual(recovery.transformForward(expected),coefficients,AbsTol=3e-12);
            end
            changedTable = sampling.observations;
            changedTable.zLower(1) = -.27; changedTable.zUpper(1) = -.27;
            replacement = IMVerticalSampling(changedTable);
            testCase.verifyNotEqual(replacement.observations,sampling.observations);
            testCase.verifyEqual(recovery.metadata,identity);
            testCase.verifyEqual(recovery.synthesisMatrix,snapshot,AbsTol=0);
            testCase.verifyEqual(counter.solveCount,initialCount);
        end

        function underflowCannotSilentlyEraseARepresentableIntegral(testCase)
            domain = [0 1e-306];
            background = IMModeBackground.constantStratification(N0=1,zDomain=domain);
            info = table("tiny:constant","intrinsic",1,VariableNames=["columnLabel","branch","branchIndex"]);
            evaluator = @(variable,z,order) double(order == 0)*(1e-16*double(variable == "F")+double(variable == "G"))*ones(numel(z),1);
            basis = IMModeBasis(background=background,family="underflowControl",modeInfo=info,normalization="maxAbsG",horizontalKind="independent",intrinsicEvaluator=evaluator,variables=["F","G"]);
            testCase.verifyGreaterThan(diff(domain)*1e-16,0);
            means = IMVerticalSampling.cellAverages("F",domain,observationLabels="tiny:mean");
            integrals = IMVerticalSampling.cellIntegrals("F",domain,observationLabels="tiny:integral");
            testCase.verifyEqual(means.evaluate(basis,initialOrder=128,maxOrder=256),1e-16,RelTol=3e-14,AbsTol=1e-30);
            testCase.verifyError(@() integrals.evaluate(basis,initialOrder=128,maxOrder=256),"IMVerticalSampling:UnrepresentableIntegral");
        end
    end
end

function observations = prescribedFixture()
root = fileparts(mfilename("fullpath"));
observations = readtable(fullfile(root,"Fixtures","IMPrescribedStaggeredSampling.csv"),TextType="string");
observations.units = string(observations.units);
end

function basis = fourierBasis(harmonics)
background = IMModeBackground.constantStratification(N0=1,zDomain=[-2 0],g=2);
columnLabel = "fourier:"+string(harmonics(:));
info = table(columnLabel,repmat("fourier",numel(harmonics),1),harmonics(:),VariableNames=["columnLabel","branch","branchIndex"]);
basis = IMModeBasis(background=background,family="samplingControl",modeInfo=info,normalization="maxAbsF",horizontalKind="independent", ...
    intrinsicEvaluator=@(variable,z,order) fourierValues(harmonics,background.zDomain,variable,z,order),variables=["F","Fz","G"]);
end

function values = fourierValues(harmonics,domain,variable,z,order)
k = pi*harmonics/diff(domain);
theta = (z(:)-domain(1))*k;
if variable == "F"
    values = cos(theta+order*pi/2).*k.^order;
else
    values = sin(theta+order*pi/2).*k.^order;
end
values(:,harmonics == 0) = 0;
if variable == "F" && order == 0, values(:,harmonics == 0) = 1; end
end

function values = exactFourierRows(observations,harmonics,domain)
k = pi*harmonics/diff(domain);
values = zeros(height(observations),numel(harmonics));
for row = 1:height(observations)
    variable = observations.variable(row); order = observations.derivativeOrder(row);
    if variable == "Fz", variable = "F"; order = order+1; end
    if variable == "Gz", variable = "G"; order = order+1; end
    if variable == "Gzz", variable = "G"; order = order+2; end
    lower = observations.zLower(row); upper = observations.zUpper(row);
    midpoint = .5*(lower+upper);
    values(row,:) = fourierValues(harmonics,domain,variable,midpoint,order);
    if observations.kind(row) ~= "point"
        width = upper-lower;
        halfPhase = .5*width*k;
        factor = ones(size(halfPhase)); nonzero = halfPhase ~= 0;
        factor(nonzero) = sin(halfPhase(nonzero))./halfPhase(nonzero);
        values(row,:) = values(row,:).*factor;
        if observations.kind(row) == "cellIntegral", values(row,:) = width*values(row,:); end
    end
end
end

function values = exactNeutralRows(basis,observations,intervals)
values = zeros(height(observations),basis.nColumns);
for column = 1:basis.nColumns
    interval = intervals(basis.modeInfo.intervalIndex(column),:);
    mode = basis.modeInfo.localModeNumber(column);
    k = pi*mode/diff(interval);
    sampled = basis.G(interval(1)+pi/(2*k));
    amplitude = sampled(column);
    for row = 1:height(observations)
        lower = max(observations.zLower(row),interval(1));
        upper = min(observations.zUpper(row),interval(2));
        if lower >= upper, continue; end
        order = observations.derivativeOrder(row); factor = 1;
        if observations.variable(row) == "APV", order = 1; factor = -basis.background.f0; end
        phaseLower = k*(lower-interval(1))+order*pi/2;
        phaseUpper = k*(upper-interval(1))+order*pi/2;
        values(row,column) = factor*amplitude*k^(order-1)*(cos(phaseLower)-cos(phaseUpper));
        if observations.kind(row) == "cellAverage"
            values(row,column) = values(row,column)/(observations.zUpper(row)-observations.zLower(row));
        end
    end
end
end

function observations = replaceWithGz(observations)
observations.variable(:) = "G";
observations.derivativeOrder(:) = 1;
end

function integral = integratePchip(knots,values,bounds)
pp = pchip(knots,values.');
[breaks,coefficients,nPieces,order,dimension] = unmkpp(pp);
integral = zeros(1,dimension);
powers = order:-1:1;
for piece = 1:nPieces
    lower = max(bounds(1),breaks(piece)); upper = min(bounds(2),breaks(piece+1));
    if lower >= upper, continue; end
    offsets = ((upper-breaks(piece)).^powers-(lower-breaks(piece)).^powers)./powers;
    rows = (piece-1)*dimension+(1:dimension);
    integral = integral+(coefficients(rows,:)*offsets.').';
end
end
