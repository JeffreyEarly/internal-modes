classdef IMModeBasisNormalizationTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function waveRulesAgreeWithFreshConstructionAcrossBackends(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            z = linspace(-1000,0,31).';
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48,coordinateKind="bandwidth")}
                original = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=3,solver=solver{1});
                originalF = original.F(z);
                for name = original.capabilities.normalizationNames
                    changed = original.withNormalization(name);
                    fresh = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=3,solver=solver{1},normalization=name);
                    testCase.verifyEqual(changed.F(z),fresh.F(z),AbsTol=3e-12,RelTol=3e-12);
                    testCase.verifyEqual(changed.G(z),fresh.G(z),AbsTol=3e-12,RelTol=3e-12);
                    testCase.verifyEqual(changed.F(z),background.D*changed.Gz(z),AbsTol=3e-13);
                    testCase.verifyEqual(changed.evaluate(z,"F",1),background.D*changed.Gzz(z),AbsTol=3e-14);
                    testCase.verifyPreserved(original,changed);
                end
                testCase.verifyEqual(original.F(z),originalF,AbsTol=0);
                testCase.verifyEqual(original.normalization,"waveEnergy");
            end
        end

        function signedEnergyRulesPreserveMetricSignatures(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            energy = IMGeneralizedEnergyBoundary(g0=-.04,gd=.03);
            solver = IMSolverSpectral(nEVP=64);
            original = IMModes.meanDensityAnomalyModes(background,energyBoundary=energy,nPositiveModes=2,solver=solver);
            testCase.verifyTrue(any(original.normSignatures < 0));
            z = linspace(-1000,0,35).';
            for name = original.capabilities.normalizationNames
                changed = original.withNormalization(name);
                fresh = IMModes.meanDensityAnomalyModes(background,energyBoundary=energy,nPositiveModes=2,solver=solver,normalization=name);
                testCase.verifyEqual(changed.G(z),fresh.G(z),AbsTol=3e-11,RelTol=3e-12);
                testCase.verifyEqual(changed.F(z),fresh.F(z),AbsTol=3e-11,RelTol=3e-12);
                restored = changed.withNormalization("generalizedEnergy");
                testCase.verifyEqual(restored.normSignatures,original.normSignatures);
                testCase.verifyEqual(restored.G(z),original.G(z),AbsTol=2e-13);
                testCase.verifyPreserved(original,changed);
            end
        end

        function remainingEVPFamiliesAgreeWithTheirFactoryRules(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            solver = IMSolverSpectral(nEVP=48);
            energy = IMGeneralizedEnergyBoundary(g0=.03,gd=.01);
            factories = {@(name) IMModes.inertialModes(background,nPositiveModes=2,solver=solver,normalization=name), ...
                @(name) IMModes.geostrophicAPVModes(background,nPositiveModes=2,solver=solver,energyBoundary=energy,normalization=name), ...
                @(name) IMModes.hydrostaticGModes(background,nInternalModes=2,solver=solver,normalization=name)};
            z = [-950;-610;-120;0];
            for factory = factories
                original = factory{1}("depthRMSF");
                for name = original.capabilities.normalizationNames
                    changed = original.withNormalization(name);
                    fresh = factory{1}(name);
                    testCase.verifyEqual(changed.F(z),fresh.F(z),AbsTol=5e-11,RelTol=5e-12);
                    testCase.verifyEqual(changed.G(z),fresh.G(z),AbsTol=5e-11,RelTol=5e-12);
                    testCase.verifyEqual(changed.Gz(z),fresh.Gz(z),AbsTol=5e-12,RelTol=5e-12);
                    testCase.verifyPreserved(original,changed);
                end
            end
        end

        function generalizedEnstrophyRulesReuseTheSameSelectedModes(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            solver = IMSolverSpectral(nEVP=48);
            original = IMModes.geostrophicGeneralizedPotentialEnstrophyModes(background,1e-3,nModes=3,solver=solver);
            selected = original.selectColumns([3 1]);
            z = [-950;-580;-100;0];
            for name = selected.capabilities.normalizationNames
                changed = selected.withNormalization(name);
                fresh = IMModes.geostrophicGeneralizedPotentialEnstrophyModes(background,1e-3,nModes=3,solver=solver,normalization=name).selectColumns([3 1]);
                testCase.verifyEqual(changed.F(z),fresh.F(z),AbsTol=3e-11,RelTol=3e-12);
                testCase.verifyEqual(changed.G(z),fresh.G(z),AbsTol=3e-11,RelTol=3e-12);
                testCase.verifyPreserved(selected,changed);
            end
        end

        function selectedHydrostaticModesDoNotNormalizeDiscardedNullColumns(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            original = IMModes.hydrostaticFModes(background,nInternalModes=3,solver=IMSolverAnalytical());
            for name = ["depthRMSG","maxAbsG","generalizedEnergy"]
                testCase.verifyError(@() original.withNormalization(name),"IMModeBasis:NullNormalization");
            end
            selected = original.selectColumns([4 2]);
            changed = selected.withNormalization("generalizedEnergy");
            testCase.verifyEqual(changed.normSignatures,[1 1]);
            z = linspace(-1000,0,401).';
            integralNorms = trapz(z,(background.N2(z)/background.g).*changed.G(z).^2);
            testCase.verifyEqual(integralNorms,[1 1],AbsTol=3e-13);
            testCase.verifyPreserved(selected,changed);
            originalF = original.F([-1000;0]);
            testCase.verifyEqual(originalF(:,1),[1;1],AbsTol=2e-15);
        end

        function neutralLabelsScaleEveryChannelAndRetainDDIdentities(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            intervals = [-1000 -750;-200 0];
            original = IMModes.neutralLabelModes(background,neutralIntervals=intervals,nModesPerInterval=[3 2]).selectColumns([5 1 4]);
            changed = original.withNormalization("maxAbsG");
            z = [-1000;-930;-800;-750;-200;-120;-10;0];
            ratio = original.scientificSpecification.resolved.normalization.scales./changed.scientificSpecification.resolved.normalization.scales;
            for variable = ["G","Gz","Gzz","APV"]
                testCase.verifyEqual(changed.evaluate(z,variable),original.evaluate(z,variable).*ratio,AbsTol=2e-14);
            end
            testCase.verifyEqual(changed.G(intervals(:)),zeros(4,3),AbsTol=0);
            testCase.verifyEqual(changed.withNormalization("depthRMSG").G(z),original.G(z),AbsTol=0);
            testCase.verifyPreserved(original,changed);
        end

        function highSelectedHarmonicsKeepTheirNormalizationResolution(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            original = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=160,solver=IMSolverAnalytical(),normalization="depthRMSF");
            columns = [160 137];
            selected = original.selectColumns(columns);
            changed = selected.withNormalization("depthRMSG");
            normalizedBeforeSelection = original.withNormalization("depthRMSG").selectColumns(columns);
            scales = changed.scientificSpecification.resolved.normalization.scales;
            expectedScales = selected.scientificSpecification.resolved.normalization.scales./(pi*columns);
            testCase.verifyEqual(scales,expectedScales,RelTol=5e-13);
            testCase.verifyEqual(scales,normalizedBeforeSelection.scientificSpecification.resolved.normalization.scales,RelTol=5e-13);
            testCase.verifyPreserved(selected,changed);
            f0 = background.f0;
            criticalBackground = IMModeBackground(N2=@(z) f0^2*ones(size(z)),zDomain=background.zDomain);
            kernel = IMModes.zeroDepthWaveModesAtWavenumber(criticalBackground,1e-3,criticalIntervals=[-900 -700],nModesPerInterval=160);
            selectedKernel = kernel.selectColumns(columns).withNormalization("depthRMSG");
            fullKernel = kernel.withNormalization("depthRMSG").selectColumns(columns);
            z = [-899;-851;-799;-741;-701];
            testCase.verifyEqual(selectedKernel.G(z),fullKernel.G(z),AbsTol=2e-12,RelTol=5e-13);
            expectedF = background.D*sqrt(2*background.D/200)*(pi*columns/200).*cos((z+900)*pi*columns/200);
            testCase.verifyEqual(selectedKernel.F(z),expectedF,AbsTol=3e-9,RelTol=5e-12);
            % Dropping the highest source column must keep the original
            % integration rule, including for source-less intrinsic modes.
            lowerColumns = [137 100];
            selectedLower = kernel.selectColumns(lowerColumns).withNormalization("depthRMSG");
            fullLower = kernel.withNormalization("depthRMSG").selectColumns(lowerColumns);
            testCase.verifyEqual(selectedLower.scientificSpecification.resolved.normalization.scales,fullLower.scientificSpecification.resolved.normalization.scales,AbsTol=0);
        end

        function fixedFrequencyRootsAndEmptyInventoriesKeepTheirContracts(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            roots = IMModes.waveModesAtFrequency(background,.5*abs(background.f0),nPropagatingModes=0,nEvanescentModes=3,solver=IMSolverAnalytical()).selectColumns([3 1]);
            changed = roots.withNormalization("maxAbsG");
            testCase.verifyLessThan(changed.kSquared,[0 0]);
            testCase.verifyEqual(changed.kSquared,roots.kSquared,AbsTol=0);
            testCase.verifyPreserved(roots,changed);
            empty = roots.selectColumns([]);
            for name = empty.capabilities.normalizationNames
                converted = empty.withNormalization(name);
                testCase.verifySize(converted.F([-1000;-500;0]),[3 0]);
                testCase.verifyEmpty(converted.normSignatures);
                testCase.verifyEqual(converted.capabilities,empty.capabilities);
                testCase.verifyPreserved(empty,converted);
            end
            neutral = IMModes.neutralLabelModes(IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]),neutralIntervals=[-1000 0],nModesPerInterval=0);
            testCase.verifySize(neutral.withNormalization("maxAbsG").APV(-500),[1 0]);
            testCase.verifyError(@() empty.withNormalization("waveEnergy"),"IMModeBasis:UnavailableNormalization");
        end

        function endpointResponsesReuseTheExistingBoundaryValueSolve(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            counter = IMGeostrophicZeroAPVTestCounter();
            solver = IMInstrumentedGeostrophicZeroAPVSolver(nEVP=48,solveCounter=counter);
            original = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1e-3,solver=solver,energyBoundary=IMGeneralizedEnergyBoundary(g0=.03,gd=.01));
            calls = counter.callCount;
            testCase.verifyGreaterThan(calls,0);
            z = [-950;-500;-20];
            for name = original.capabilities.normalizationNames
                changed = original.withNormalization(name);
                ratio = original.scientificSpecification.resolved.normalization.scales./changed.scientificSpecification.resolved.normalization.scales;
                testCase.verifyEqual(changed.F(z),original.F(z).*ratio,AbsTol=3e-13);
                testCase.verifyEqual(changed.G(z),original.G(z).*ratio,AbsTol=3e-13);
                testCase.verifyPreserved(original,changed);
            end
            testCase.verifyEqual(counter.callCount,calls);
        end

        function sameNameAndInvalidRequestsAvoidProfileEvaluation(testCase)
            calls = 0;
            background = IMModeBackground(N2=@(z) ones(size(z)),zDomain=[-1 0]);
            info = table("intrinsic:1","intrinsic",1,VariableNames=["columnLabel","branch","branchIndex"]);
            original = IMModeBasis(background=background,family="intrinsic",modeInfo=info,normalization="depthRMSF",intrinsicEvaluator=@evaluate,capabilities="F");
            calls = 0;
            changed = original.withNormalization("depthRMSF");
            testCase.verifyEqual(changed.scientificSpecification,original.scientificSpecification);
            testCase.verifyEqual(calls,0);
            testCase.verifyError(@() original.withNormalization("unknownRule"),"IMModeBasis:UnknownNormalization");
            testCase.verifyError(@() original.withNormalization("depthRMSG"),"IMModeBasis:UnavailableNormalization");
            testCase.verifyEqual(calls,0);

            function values = evaluate(~,z,~)
                calls = calls+1;
                values = ones(size(z));
            end
        end

        function nonfiniteMagnitudeIsRejectedWithoutChangingTheInput(testCase)
            background = IMModeBackground(N2=@(z) ones(size(z)),zDomain=[-1 0]);
            info = table("intrinsic:1","intrinsic",1,VariableNames=["columnLabel","branch","branchIndex"]);
            original = IMModeBasis(background=background,family="intrinsic",modeInfo=info,normalization="depthRMSG",intrinsicEvaluator=@evaluate,variables=["F","G"]);
            testCase.verifyError(@() original.withNormalization("depthRMSF"),"IMModeBasis:NullNormalization");
            testCase.verifyEqual(original.G(-.5),1,AbsTol=2e-14);
            testCase.verifyEqual(original.normalization,"depthRMSG");

            function values = evaluate(variable,z,~)
                if variable == "F", values = 1e308*ones(size(z)); else, values = ones(size(z)); end
            end
        end

        function normalizationDoesNotEnterEVPOrWeightFittingConstruction(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            original = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=3,solver=IMSolverSpectral(nEVP=48));
            cleanup = onCleanup(@() profile("off"));
            profile("clear");
            profile("on");
            changed = original.withNormalization("depthRMSG");
            profile("off");
            evidence = profile("info");
            names = string({evidence.FunctionTable.FunctionName});
            forbidden = ["solveRequestedEVP","solveAssembledEVP","solveEVP","solveGeostrophicZeroAPVModes","quadratureWeightsForPoints","IMQuadratureWeightFit","fitDiscreteTransform"];
            testCase.verifyFalse(any(contains(names,forbidden)));
            testCase.verifyTrue(any(contains(names,"withNormalization")));
            testCase.verifyPreserved(original,changed);
            clear cleanup
        end
    end

    methods (Access = private)
        function verifyPreserved(testCase,original,changed)
            testCase.verifyEqual(changed.family,original.family);
            testCase.verifyEqual(changed.columnLabels,original.columnLabels);
            testCase.verifyEqual(changed.modeInfo,original.modeInfo);
            testCase.verifyEqual(changed.horizontalKind,original.horizontalKind);
            testCase.verifyEqual(changed.k,original.k);
            testCase.verifyEqual(changed.kSquared,original.kSquared);
            testCase.verifyEqual(changed.frequency,original.frequency);
            testCase.verifyEqual(changed.capabilities,original.capabilities);
            testCase.verifyEqual(changed.metadata,original.metadata);
            testCase.verifyEqual(changed.constructionRecord,original.constructionRecord);
            expected = original.scientificSpecification;
            expected.resolved.normalization = changed.scientificSpecification.resolved.normalization;
            if original.family == "geostrophicZeroAPV"
                ratio = original.scientificSpecification.resolved.normalization.scales./changed.scientificSpecification.resolved.normalization.scales;
                expectedMap = original.sourceCoordinateMap.*ratio;
                testCase.verifyEqual(changed.sourceCoordinateMap,expectedMap,AbsTol=3e-14);
                expected.resolved.sourceCoordinateMap = expectedMap;
            end
            testCase.verifyEqual(changed.scientificSpecification,expected);
            testCase.verifyEqual(changed.scientificSpecification.resolved.normalization.name,changed.normalization);
            testCase.verifyEqual(changed.scientificSpecification.resolved.normalization.normSignatures,changed.normSignatures);
            for name = ["h","eigenvalues"]
                if ismember(name,original.capabilities.spectrumNames)
                    testCase.verifyEqual(changed.(name),original.(name));
                end
            end
        end
    end
end
