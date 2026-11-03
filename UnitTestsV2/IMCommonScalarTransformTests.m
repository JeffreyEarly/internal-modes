classdef IMCommonScalarTransformTests < matlab.unittest.TestCase
    % Check physical source pairings independently of sample recovery.

    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function apvPairingsMatchIndependentFormsAndNonmodalSources(testCase)
            background = signedBackground();
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            solvers = {IMSolverAnalytical(),IMSolverSpectral(nEVP=64,coordinateKind="z")};
            [z,w] = endpointRule(background.zDomain,96);
            for j = 1:numel(solvers)
                basis = IMModes.geostrophicAPVModes(background,energyBoundary=energy,nPositiveModes=3,solver=solvers{j});
                for variable = ["F","G"]
                    physical = basis.projectionOnGrid(z,w,variable=variable,prepareRecovery=false);
                    [H,Hplus] = independentForms(basis,variable,energy,192);
                    [Href,~] = independentForms(basis,variable,energy,256);
                    testCase.verifyEqual(H,Href,RelTol=3e-12,AbsTol=3e-12);
                    tolerance = 3e-9;
                    if j == 1, tolerance = 3e-5; end
                    C = chol(Hplus)\eye(basis.nColumns);
                    targetError = norm(C'*(physical.targetGramMatrix-H)*C,2);
                    majorantError = norm(C'*(physical.majorantGramMatrix-Hplus)*C,2);
                    fprintf("Common APV %s/%s physical operator errors: H %.3g, Hplus %.3g.\n",string(class(solvers{j})),variable,targetError,majorantError);
                    testCase.verifyLessThanOrEqual(targetError,tolerance);
                    testCase.verifyLessThanOrEqual(majorantError,tolerance);
                    A = basis.evaluate(z,variable);
                    metric = independentMetric(background,energy,"geostrophicAPV",variable,z,w);
                    P = A'*diag(metric);
                    testCase.verifyEqual(physical.samplePairingMatrix,P,AbsTol=2e-13);
                    testCase.verifyEqual(physical.columnLabels,basis.columnLabels);
                    testCase.verifyTrue(all(physical.activeColumnMask));
                    input = independentSource(z,variable);
                    expected = physical.targetGramMatrix\(P*input);
                    testCase.verifyEqual(physical.projectPhysical(input),expected,RelTol=2e-12,AbsTol=2e-12);
                    [zr,wr] = endpointRule(background.zDomain,256);
                    refinedMetric = independentMetric(background,energy,"geostrophicAPV",variable,zr,wr);
                    continuous = Href\(basis.evaluate(zr,variable)'*(refinedMetric.*independentSource(zr,variable)));
                    testCase.verifyEqual(physical.projectPhysical(input),continuous,RelTol=2*tolerance,AbsTol=2*tolerance);
                    testCase.verifyEmpty(physical.recovery);
                    if variable == "G"
                        testCase.verifyTrue(any(eig(H) < -.01));
                        testCase.verifyGreaterThan(min(eig(Hplus)),0);
                    end
                end
            end
        end

        function literalEndpointCoefficientsSurvivePhysicalUnitScaling(testCase)
            [z,w] = endpointRule([-1 0],16);
            for g = [2 1e-16]
                background = IMModeBackground.constantStratification(N0=sqrt(g),zDomain=[-1 0],g=g,rotationRate=.5,latitude=90,surfaceCondition="freeSurface");
                for g0 = [.3 1e-15]
                    boundary = IMGeneralizedEnergyBoundary(g0=g0,gd=Inf);
                    basis = IMModes.geostrophicAPVModes(background,energyBoundary=boundary,nPositiveModes=0,solver=IMSolverAnalytical());
                    physical = basis.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
                    expected = g0/(g+g0);
                    testCase.verifyEqual(physical.metricMatrix(end,end),expected,RelTol=1e-14);
                    testCase.verifyError(@() basis.projectionOnGrid(z(1:end-1),w(1:end-1),variable="G",prepareRecovery=false),"IMProjection:MissingEndpointSample");
                    mda = IMModes.meanDensityAnomalyModes(background,energyBoundary=boundary,nPositiveModes=0,solver=IMSolverSpectral(nEVP=32,coordinateKind="z"));
                    physical = mda.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
                    testCase.verifyEqual(physical.metricMatrix(end,end),g0/g,RelTol=1e-14);
                end
            end
        end

        function mdaUsesItsOwnSignedGFormAndExplicitRecoveryPolicy(testCase)
            background = signedBackground();
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            basis = IMModes.meanDensityAnomalyModes(background,energyBoundary=energy,nPositiveModes=3,solver=IMSolverSpectral(nEVP=64,coordinateKind="z"));
            [z,w] = endpointRule(background.zDomain,96);
            [H,Hplus] = independentForms(basis,"G",energy,256);
            transform = basis.scalarTransform(z,w,variable="G");
            testCase.verifyClass(transform,"IMScalarDiscreteTransform");
            testCase.verifyTrue(transform.success);
            testCase.verifyEqual(transform.targetGramMatrix,H,RelTol=3e-9,AbsTol=3e-9);
            testCase.verifyEqual(transform.projection.majorantGramMatrix,Hplus,RelTol=3e-9,AbsTol=3e-9);
            testCase.verifyTrue(any(eig(H) < -.01));
            testCase.verifyEqual(transform.columnLabels,basis.columnLabels);
            testCase.verifyEmpty(transform.modeNumber);
            testCase.verifyEqual(transform.z,z,AbsTol=0);
            testCase.verifyEqual(transform.weights,w,AbsTol=0);
            A = basis.G(z);
            metric = independentMetric(background,energy,"meanDensityAnomaly","G",z,w);
            B = transform.targetGramMatrix\(A'*diag(metric));
            source = independentSource(z,"G");
            testCase.verifyEqual(transform.physicalDualMatrix,B,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(transform.projection.projectPhysical(source),B*source,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(transform.transformForward(A),eye(basis.nColumns),AbsTol=3e-11);
            testCase.verifyError(@() basis.scalarTransform(z,w,variable="F"),"IMModeBasis:UnavailableProjection");
            testCase.verifyError(@() basis.scalarTransform(z,w,variable="G",policy="sampledLeastSquares"),"IMModeBasis:MissingFitFactor");
            fit = diag(linspace(.6,1.4,numel(z)));
            leastSquares = basis.scalarTransform(z,w,variable="G",policy="sampledLeastSquares",fitFactor=fit);
            testCase.verifyTrue(leastSquares.success);
            testCase.verifyEqual(leastSquares.forwardMatrix,(fit*A)\fit,RelTol=2e-12,AbsTol=2e-12);
            testCase.verifyEqual(leastSquares.physicalDualMatrix,B,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyGreaterThan(norm((leastSquares.forwardMatrix-transform.forwardMatrix)*source),1e-5);
        end

        function requiredEndpointTracesAndInactiveLimitsAreLiteral(testCase)
            background = signedBackground();
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            basis = IMModes.geostrophicAPVModes(background,energyBoundary=energy,nPositiveModes=2,solver=IMSolverSpectral(nEVP=48,coordinateKind="z"));
            [z,w] = endpointRule(background.zDomain,48);
            testCase.verifyError(@() basis.scalarTransform(z(1:end-1),w(1:end-1),variable="G"),"IMProjection:MissingEndpointSample");
            testCase.verifyError(@() basis.scalarTransform(z(2:end),w(2:end),variable="G"),"IMProjection:MissingEndpointSample");
            interior = basis.scalarTransform(z(2:end-1),w(2:end-1),variable="F");
            testCase.verifyTrue(interior.success);
            rigid = signedBackground("rigidLid");
            noTraces = IMModes.geostrophicAPVModes(rigid,energyBoundary=IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf),nPositiveModes=2,solver=IMSolverSpectral(nEVP=48,coordinateKind="z"));
            unaugmented = noTraces.scalarTransform(z(2:end-1),w(2:end-1),variable="G");
            testCase.verifyTrue(unaugmented.success);
            testCase.verifyEqual(unaugmented.metricMatrix,diag(w(2:end-1).*rigid.N2(z(2:end-1))/rigid.g),AbsTol=0);
            freeInfinite = IMModes.geostrophicAPVModes(background,energyBoundary=IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf),nPositiveModes=2,solver=IMSolverSpectral(nEVP=48,coordinateKind="z"));
            testCase.verifyError(@() freeInfinite.scalarTransform(z(2:end-1),w(2:end-1),variable="G"),"IMProjection:MissingEndpointSample");
            freePhysical = freeInfinite.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            testCase.verifyEqual(freePhysical.metricMatrix(end,end),1,AbsTol=0);
            singular = IMModes.geostrophicAPVModes(background,energyBoundary=IMGeneralizedEnergyBoundary(g0=-background.g,gd=Inf),nPositiveModes=2,solver=IMSolverSpectral(nEVP=48,coordinateKind="z"));
            singularPhysical = singular.projectionOnGrid(z(2:end-1),w(2:end-1),variable="G",prepareRecovery=false);
            testCase.verifyEqual(singularPhysical.metricMatrix,diag(w(2:end-1).*background.N2(z(2:end-1))/background.g),AbsTol=0);
        end

        function unsupportedDerivativeAndCompanionObservationsStayUnsupported(testCase)
            background = signedBackground();
            solver = IMSolverSpectral(nEVP=48,coordinateKind="z");
            evp = IMInternalModes(zDomain=background.zDomain,N2=@(z) ones(size(z)),formulation="G",modeFamily="hydrostatic",r=.5,rotationRate=background.rotationRate,latitude=background.latitude,g=background.g, ...
                surfaceBoundary=IMBoundaryCondition(a=1,b=0,c=0,d=1),bottomBoundary=IMBoundaryCondition.dirichlet());
            source = solver.solveEVP(evp,nModes=3);
            info = table("derivative:"+string((1:3).'),repmat("positive",3,1),(1:3).',VariableNames=["columnLabel","branch","branchIndex"]);
            wrapper = IMModeBasis(source=source,background=background,family="hydrostaticG",modeInfo=info,normalization="depthRMSG",horizontalKind="independent");
            [z,w] = endpointRule(background.zDomain,48);
            testCase.verifySize(wrapper.G(z),[numel(z) 3]);
            testCase.verifyError(@() wrapper.scalarTransform(z,w,variable="G"),"IMModeBasis:UnavailableProjection");
            testCase.verifyError(@() wrapper.selectColumns([]).scalarTransform(z,w,variable="G"),"IMModeBasis:UnavailableProjection");
            coupled = IMModes.geostrophicGeneralizedPotentialEnstrophyModes(background,.8,nModes=2,alpha0=.7,alphaD=.4,solver=solver);
            testCase.verifySize(coupled.G(z),[numel(z) coupled.nColumns]);
            testCase.verifyError(@() coupled.scalarTransform(z,w,variable="G"),"IMModeBasis:UnavailableProjection");
            waveBackground = IMModeBackground.constantStratification(N0=2,zDomain=[-1 0],g=2,rotationRate=background.rotationRate,latitude=33);
            wave = IMModes.waveModesAtWavenumber(waveBackground,.3,nInternalModes=2,solver=solver);
            testCase.verifyError(@() wave.scalarTransform(z,w,variable="F"),"IMModeBasis:UnavailableProjection");
            zero = IMModes.geostrophicZeroAPVModesAtWavenumber(background,.3,energyBoundary=IMGeneralizedEnergyBoundary(g0=0,gd=0),solver=solver);
            testCase.verifyError(@() zero.scalarTransform(z,w,variable="G"),"IMModeBasis:UnavailableProjection");
        end

        function selectedNormalizationCovariesWithFixedRecoveryFactors(testCase)
            background = signedBackground();
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            full = IMModes.geostrophicAPVModes(background,energyBoundary=energy,nPositiveModes=3,solver=IMSolverSpectral(nEVP=64,coordinateKind="z"));
            columns = [full.nColumns 1 2];
            basis = full.selectColumns(columns);
            changed = basis.withNormalization("maxAbsG");
            T = normalizationMap(basis,changed);
            [z,w] = endpointRule(background.zDomain,64);
            physical = basis.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            changedPhysical = changed.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            fullPhysical = full.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            testCase.verifyEqual(physical.targetGramMatrix,fullPhysical.targetGramMatrix(columns,columns),AbsTol=0);
            testCase.verifyEqual(physical.samplePairingMatrix,fullPhysical.samplePairingMatrix(columns,:),AbsTol=0);
            testCase.verifyEqual(physical.physicalDualMatrix,physical.targetGramMatrix\physical.samplePairingMatrix,RelTol=3e-12,AbsTol=3e-12);
            C = chol(physical.majorantGramMatrix)\eye(basis.nColumns);
            D = diag(linspace(.8,1.3,numel(z)));
            F = diag(linspace(1.4,.6,numel(z)));
            original = basis.scalarTransform(z,w,variable="G",coefficientScale=C,observationFactor=D,fitFactor=F);
            normalized = changed.scalarTransform(z,w,variable="G",coefficientScale=T\C,observationFactor=D,fitFactor=F);
            testCase.verifyTrue(original.success && normalized.success);
            testCase.verifyEqual(changed.G(z),basis.G(z)*T,RelTol=3e-14,AbsTol=3e-14);
            testCase.verifyEqual(changedPhysical.targetGramMatrix,T'*physical.targetGramMatrix*T,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(changedPhysical.majorantGramMatrix,T'*physical.majorantGramMatrix*T,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(changedPhysical.samplePairingMatrix,T'*physical.samplePairingMatrix,RelTol=3e-14,AbsTol=3e-14);
            testCase.verifyEqual(normalized.physicalDualMatrix,T\original.physicalDualMatrix,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(normalized.forwardMatrix,T\original.forwardMatrix,RelTol=3e-12,AbsTol=3e-12);
            input = independentSource(z,"G");
            testCase.verifyEqual(normalized.transformForward(input),T\original.transformForward(input),RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(normalized.transformBack(normalized.transformForward(input)),original.transformBack(original.transformForward(input)),RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(normalized.recoverySpecification.observationFactor,D,AbsTol=0);
            testCase.verifyEqual(normalized.recoverySpecification.fitFactor,F,AbsTol=0);
            testCase.verifyEqual(normalized.recoverySpecification.coefficientScale,T\C,AbsTol=0);
            testCase.verifyEqual(original.columnLabels,full.columnLabels(columns));
            supplied = IMDiscreteTransform.fromObservationMatrices(basis.G(z),policy="sampledLeastSquares",columnLabels=basis.columnLabels,fitFactor=F,physicalDual=physical.physicalDualMatrix);
            rebound = basis.scalarTransform(z,w,variable="G",recovery=supplied);
            testCase.verifyEqual(rebound.forwardMatrix,supplied.forwardMatrix,AbsTol=0);
            testCase.verifyEqual(rebound.recoverySpecification,supplied.recoverySpecification);
            collection = IMBasisCollection({basis,changed},basisIndex=[2 1 2]);
            for page = [3 2 1]
                expected = collection.bases{collection.basisIndex(page)}.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
                actual = collection.projectionOnGrid(z,w,page=page,variable="G",prepareRecovery=false);
                testCase.verifyEqual(actual.targetGramMatrix,expected.targetGramMatrix,AbsTol=0);
                testCase.verifyEqual(actual.physicalDualMatrix,expected.physicalDualMatrix,AbsTol=0);
                testCase.verifyEqual(actual.columnLabels,expected.columnLabels);
            end
        end

        function hydrostaticActivityEmptyAndFailedRecoveryRemainScientific(testCase)
            background = signedBackground("rigidLid");
            full = IMModes.hydrostaticFModes(background,nInternalModes=3,solver=IMSolverAnalytical());
            [z,w] = endpointRule(background.zDomain,64);
            for variable = ["F","G"]
                transform = full.scalarTransform(z,w,variable=variable);
                testCase.verifyTrue(transform.success);
                active = true(1,full.nColumns);
                if variable == "G", active(full.modeInfo.branch == "barotropic") = false; end
                testCase.verifyEqual(transform.activeColumnMask,active);
                testCase.verifyEqual(transform.forwardMatrix*transform.inverseMatrix,diag(double(active)),AbsTol=3e-12);
                empty = full.selectColumns([]).scalarTransform(z,w,variable=variable);
                explicitEmpty = full.projectionOnGrid(z,w,variable=variable,columns=[],prepareRecovery=false);
                testCase.verifySize(empty.forwardMatrix,[0 numel(z)]);
                testCase.verifySize(empty.inverseMatrix,[numel(z) 0]);
                testCase.verifySize(empty.targetGramMatrix,[0 0]);
                testCase.verifyEmpty(empty.columnLabels);
                testCase.verifySize(explicitEmpty.physicalDualMatrix,[0 numel(z)]);
            end
            endpoints = background.zDomain.';
            failed = full.scalarTransform(endpoints,[.5;.5],variable="G");
            testCase.verifyFalse(failed.success);
            testCase.verifyEqual(failed.activeColumnMask,[false true true true]);
            testCase.verifyEqual(failed.projection.physicalDualStatus,"available");
            testCase.verifySize(failed.physicalDualMatrix,[full.nColumns 2]);
            input = [1+2i;-.4+.7i];
            testCase.verifyEqual(failed.projection.projectPhysical(input),failed.physicalDualMatrix*input,AbsTol=0);
            testCase.verifyEmpty(failed.forwardMatrix);
            mda = IMModes.meanDensityAnomalyModes(background,energyBoundary=IMGeneralizedEnergyBoundary(g0=0,gd=0),nPositiveModes=2,solver=IMSolverSpectral(nEVP=48,coordinateKind="z"));
            constant = find(isinf(mda.h));
            testCase.verifyEqual(numel(constant),1);
            mdaTransform = mda.scalarTransform(z,w,variable="G");
            testCase.verifyTrue(mdaTransform.success);
            testCase.verifyTrue(mdaTransform.activeColumnMask(constant));
            testCase.verifyGreaterThan(norm(mdaTransform.inverseMatrix(:,constant)),1);
            free = signedBackground();
            for factory = {@IMModes.hydrostaticFModes,@IMModes.hydrostaticGModes}
                basis = factory{1}(free,nInternalModes=2,solver=IMSolverAnalytical());
                for variable = ["F","G"]
                    transform = basis.scalarTransform(z,w,variable=variable);
                    testCase.verifyTrue(transform.success);
                    testCase.verifyEqual(transform.columnLabels,basis.columnLabels);
                    testCase.verifyEqual(transform.transformForward(basis.evaluate(z,variable)),eye(basis.nColumns),AbsTol=3e-11);
                end
            end
        end

        function neutralDualAndFixedRuleReuseDoNotResolveTheEVP(testCase)
            background = signedBackground();
            counter = IMScalarTransformTestCounter();
            solver = IMScalarTransformCountingSolver(counter);
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            full = IMModes.geostrophicAPVModes(background,energyBoundary=energy,nPositiveModes=3,solver=solver);
            initialCount = counter.solveCount;
            testCase.verifyGreaterThan(initialCount,0);
            [z,w] = endpointRule(background.zDomain,64);
            selected = full.selectColumns([3 1]).withNormalization("maxAbsG");
            projection = selected.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            transform = selected.scalarTransform(z,w,variable="G");
            testCase.verifyTrue(transform.success);
            testCase.verifyEqual(counter.solveCount,initialCount);
            testCase.verifyEqual(transform.z,z,AbsTol=0);
            testCase.verifyEqual(transform.weights,w,AbsTol=0);
            testCase.verifyEqual(transform.physicalDualMatrix,projection.physicalDualMatrix,AbsTol=0);
            neutralBackground = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1 0],rotationRate=background.rotationRate,latitude=33,g=2);
            neutral = IMModes.neutralLabelModes(neutralBackground,neutralIntervals=[-1 0],energyBoundary=IMGeneralizedEnergyBoundary(g0=0,gd=0),nModesPerInterval=3);
            physical = neutral.projectionOnGrid(z,w,variable="APV",prepareRecovery=false);
            adapter = neutral.scalarTransform(z,w,variable="APV");
            testCase.verifyTrue(adapter.success);
            testCase.verifyEqual(adapter.physicalDualMatrix,physical.physicalDualMatrix,AbsTol=0);
            testCase.verifyEqual(adapter.projection.projectionKind,"prescribedDual");
            testCase.verifyFalse(adapter.projection.supportsGramAssessment);
            testCase.verifyEqual(adapter.inverseMatrix,neutral.APV(z),AbsTol=0);
            testCase.verifyEqual(adapter.transformForward(neutral.APV(z)),eye(neutral.nColumns),AbsTol=3e-12);
            input = independentSource(z,"F");
            testCase.verifyEqual(adapter.projection.projectPhysical(input),physical.physicalDualMatrix*input,AbsTol=0);
            empty = neutral.selectColumns([]).scalarTransform(z,w,variable="APV");
            testCase.verifySize(empty.forwardMatrix,[0 numel(z)]);
            testCase.verifySize(empty.physicalDualMatrix,[0 numel(z)]);
        end
    end
end

function background = signedBackground(surfaceCondition)
arguments
    surfaceCondition (1,1) string = "freeSurface"
end
background = IMModeBackground.constantStratification(N0=1,zDomain=[-1 0],g=2,rotationRate=1/(2*sind(33)),latitude=33,surfaceCondition=surfaceCondition);
end

function [z,w] = endpointRule(domain,n)
[interior,weights] = legpts(n,domain);
z = [domain(1);interior;domain(2)];
w = [0;weights(:);0];
end

function values = independentSource(z,variable)
if variable == "F"
    values = [exp(z)+.2*z.^3+1i*cos(2.3*z),sin(1.7*z)+1i*(1+z.^2)];
else
    values = [cos(1.9*z)+.4*z+1i*exp(.7*z),1+z.^3+1i*sin(2.7*z)];
end
end

function [H,Hplus] = independentForms(basis,variable,energy,n)
[z,w] = endpointRule(basis.zDomain,n);
A = basis.evaluate(z,variable);
metric = independentMetric(basis.background,energy,basis.family,variable,z,w);
H = A'*(metric.*A);
Hplus = A'*(abs(metric).*A);
end

function metric = independentMetric(background,energy,family,variable,z,w)
if variable == "F"
    metric = w;
    return
end
metric = w.*background.N2(z)/background.g;
if family == "meanDensityAnomaly"
    surface = energy.g0/background.g;
elseif background.surfaceCondition == "freeSurface"
    if isinf(energy.g0)
        surface = 1;
    elseif energy.g0 == -background.g
        surface = 0;
    else
        surface = energy.g0/(background.g+energy.g0);
    end
else
    surface = energy.g0/background.g;
end
bottom = energy.gd/background.g;
if isfinite(surface), metric(end) = metric(end)+surface; end
if isfinite(bottom), metric(1) = metric(1)+bottom; end
end

function T = normalizationMap(basis,changed)
old = basis.scientificSpecification.resolved.normalization.scales;
new = changed.scientificSpecification.resolved.normalization.scales;
T = diag(old./new);
end
