classdef IMBasisCapabilityTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function waveDerivativesAgreeAcrossBackendsAndRepeatedPages(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48,coordinateKind="z")}
                basis = IMModes.waveModesAtWavenumber(background,1e-3,solver=solver{1},nInternalModes=3);
                collection = IMBasisCollection({basis},basisIndex=[1 1 1]);
                testCase.verifyEqual(collection.metadata.capabilities,basis.capabilities);
                testCase.verifyEqual(basis.capabilities.variables,["F","G"]);
                testCase.verifyEqual(basis.capabilities.derivativeOrders,{0:1,0:2});
                z = linspace(-1000,0,21).';
                pages = [3 1 3 2]; columns = [3 1 3];
                Fz = collection.evaluate(z,variable="F",derivativeOrder=1,pages=pages,columns=columns,sampleChunkSize=4,pageChunkSize=2);
                Gzz = collection.evaluate(z,variable="Gz",derivativeOrder=1,pages=pages,columns=columns,sampleChunkSize=3,pageChunkSize=1);
                testCase.verifyEqual(Fz,background.D*Gzz,AbsTol=3e-14);
                scalar = basis.evaluate(z,"F",1);
                testCase.verifyEqual(Fz,repmat(scalar(:,columns),1,1,numel(pages)),AbsTol=3e-17);
                testCase.verifyEqual(collection.evaluate(z,variable="Gzz"),collection.evaluate(z,variable="G",derivativeOrder=2),AbsTol=0);
                testCase.verifyEqual(collection.evaluate(z,variable="Fz"),collection.evaluate(z,variable="F",derivativeOrder=1),AbsTol=0);
            end
        end

        function projectionAndRecoveryCapabilitiesRemainSeparate(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            hydrostatic = IMModes.hydrostaticFModes(background,solver=IMSolverAnalytical(),nInternalModes=2);
            wave = IMModes.waveModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),nInternalModes=2);
            testCase.verifyEqual([hydrostatic.capabilities.projectionOperations.variable],["F","G"]);
            testCase.verifyEqual([hydrostatic.capabilities.projectionOperations.kind],["sampledGalerkin","sampledGalerkin"]);
            testCase.verifyEmpty(wave.capabilities.projectionOperations);
            for basis = {hydrostatic,wave}
                description = basis{1}.capabilities;
                testCase.verifyEqual(description.sampledRecovery.availability,"requiresPreparedSystem");
                testCase.verifyEqual(description.sampledRecovery.noiseEvidence,"requiresDeclaredObservationFactor");
                testCase.verifyEqual(description.sampledRecovery.synthesisVariables,description.variables);
                testCase.verifyFalse(isfield(description.sampledRecovery,"rank"));
                testCase.verifyFalse(isfield(description.sampledRecovery,"physicalGramRank"));
            end
            % Profiles can be sampled without supplying enough information
            % to recover all retained coefficients.
            sampled = wave.F(-500);
            transform = IMDiscreteTransform.fromObservationMatrices(sampled,columnLabels=wave.columnLabels,observationFactor=1,failureMode="diagnostic");
            testCase.verifyFalse(transform.success);
            testCase.verifyEmpty(transform.forwardMatrix);
        end

        function neutralTracesAndAPVDualRequirementsAreExplicit(testCase)
            background = IMModeBackground(N2=@(z) 1e-4*double(z < -200),zDomain=[-1000 0]);
            basis = IMModes.neutralLabelModes(background,neutralIntervals=[-200 0],nModesPerInterval=2);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            description = basis.capabilities;
            testCase.verifyEqual(description.variables,["G","APV"]);
            testCase.verifyEqual(description.derivativeOrders,{0:2,0});
            testCase.verifyEqual(description.traceOrders,description.derivativeOrders);
            testCase.verifyEqual(description.derivativeConvention,"classicalIntervalInterior");
            testCase.verifyEqual(description.traceConvention,"insideIntervalEndpoints");
            testCase.verifyEmpty(description.spectrumNames(ismember(description.spectrumNames,["h","eigenvalues"])));
            testCase.verifyEqual(description.projectionOperations.variable,"APV");
            testCase.verifyEqual(description.projectionOperations.kind,"prescribedContinuousDual");
            testCase.verifyEqual(description.projectionOperations.availability,"conditional");
            testCase.verifyTrue(any(contains(description.projectionOperations.requirements,"nonzero f0")));
            testCase.verifyEqual(collection.metadata.capabilities,description);
            z = [-201;-200;-100;0];
            derivative = basis.Gz(z);
            testCase.verifyEqual(derivative(1,:),[0 0]);
            testCase.verifyGreaterThan(abs(derivative(2,:)),[0 0]);
            testCase.verifyEqual(collection.evaluate(z,variable="Gz",pages=[2 1 2]),repmat(derivative,1,1,3),AbsTol=0);
            testCase.verifyError(@() collection.evaluate(z,variable="APV",derivativeOrder=1),"IMBasisCollection:UnsupportedDerivative");
            testCase.verifyError(@() basis.evaluate(z,"APV",1),"IMModeBasis:UnsupportedDerivative");
        end

        function emptySelectionsKeepTheirDeclaredOperations(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            full = IMModes.waveModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),nInternalModes=2);
            empty = full.selectColumns([]);
            collection = IMBasisCollection({empty},basisIndex=[1 1]);
            testCase.verifyEqual(empty.capabilities,full.capabilities);
            testCase.verifyEqual(empty.capabilities.spectrumNames,["h","eigenvalues","frequency"]);
            testCase.verifyTrue(empty.capabilities.supportsEmptyColumns);
            testCase.verifySize(collection.evaluate([-1000;-500;0],variable="F",derivativeOrder=1),[3 0 2]);
            testCase.verifyError(@() empty.evaluate(-500,"F",2),"IMModeBasis:UnsupportedDerivative");
            testCase.verifyError(@() collection.evaluate(-500,variable="F",derivativeOrder=2),"IMBasisCollection:UnsupportedDerivative");
            neutralBackground = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            neutral = IMModes.neutralLabelModes(neutralBackground,neutralIntervals=[-1000 0],nModesPerInterval=0);
            testCase.verifyEmpty(neutral.capabilities.spectrumNames(ismember(neutral.capabilities.spectrumNames,["h","eigenvalues"])));
            testCase.verifyEqual(neutral.capabilities.projectionOperations.variable,"APV");
            testCase.verifyError(@() neutral.evaluate(-500,"G",3),"IMModeBasis:UnsupportedDerivative");
        end

        function unsupportedRequestsDoNotEvaluateIntrinsicProfiles(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            calls = 0;
            info = table("profile:1","intrinsic",1,VariableNames=["columnLabel","branch","branchIndex"]);
            basis = IMModeBasis(background=background,family="intrinsic",modeInfo=info,normalization="depthRMSF",intrinsicEvaluator=@evaluate,capabilities="F");
            calls = 0;
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            testCase.verifyEqual(calls,0);
            testCase.verifyEqual(basis.capabilities.derivativeOrders,{0});
            testCase.verifyError(@() basis.evaluate(-500,"F",1),"IMModeBasis:UnsupportedDerivative");
            testCase.verifyError(@() collection.evaluate(-500,variable="F",derivativeOrder=1),"IMBasisCollection:UnsupportedDerivative");
            testCase.verifyError(@() collection.productFactor(id="bad",variable="F",derivativeOrder=1),"IMBasisCollection:UnsupportedFactor");
            testCase.verifyEqual(calls,0);

            function values = evaluate(~,z,~)
                calls = calls+1;
                values = ones(size(z));
            end
        end

        function productFactorsPreserveTheSameDerivativeContract(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            basis = IMModes.waveModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),nInternalModes=2);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            Fz = collection.productFactor(id="Fz",page=2,variable="F",derivativeOrder=1);
            Gzz = collection.productFactor(id="Gzz",page=2,variable="Gz",derivativeOrder=1);
            z = [-900;-600;-100];
            columns = [3 1];
            testCase.verifyEqual(Fz.evaluate(z,columns,"reference"),background.D*Gzz.evaluate(z,columns,"reference"),AbsTol=2e-14);
            testCase.verifyEqual(Gzz.provenance.variable,"G");
            testCase.verifyEqual(Gzz.provenance.derivativeOrder,2);
            testCase.verifyError(@() collection.productFactor(id="bad",variable="Gzz",derivativeOrder=1),"IMBasisCollection:UnsupportedFactor");
        end

        function fixedFrequencyDoesNotClaimAFixedWavenumberEnergyNorm(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            basis = IMModes.waveModesAtFrequency(background,.005,solver=IMSolverAnalytical(),nPropagatingModes=2);
            testCase.verifyTrue(ismember(basis.normalization,basis.capabilities.normalizationNames));
            testCase.verifyFalse(ismember("waveEnergy",basis.capabilities.normalizationNames));
            testCase.verifyEqual(basis.capabilities.spectrumNames,["h","eigenvalues","kSquared","frequency"]);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            testCase.verifyEqual(collection.metadata.capabilities,basis.capabilities);
            z = [-900;-100];
            testCase.verifyEqual(collection.evaluate(z,variable="F",derivativeOrder=1),repmat(basis.evaluate(z,"F",1),1,1,2),AbsTol=0);
        end

        function endpointResponsesExposeTheirRegularPairDerivatives(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            energy = IMGeneralizedEnergyBoundary(g0=0,gd=0);
            basis = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),energyBoundary=energy);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            z = [-1000;-600;0];
            Fz = collection.evaluate(z,variable="F",derivativeOrder=1);
            G = collection.evaluate(z,variable="G");
            testCase.verifyEqual(Fz,-background.N2(z)/background.g.*G,AbsTol=2e-15);
            Gz = collection.evaluate(z,variable="Gz");
            F = collection.evaluate(z,variable="F");
            testCase.verifyEqual(Gz,-background.g*1e-6/background.f0^2.*F,AbsTol=2e-13);
            testCase.verifyEmpty(basis.capabilities.projectionOperations);
            testCase.verifyTrue(ismember("endpointResponse",basis.capabilities.normalizationNames));
        end

        function requestedEmptyAndSelectedEmptyPreserveEVPNormalizations(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            solver = IMSolverSpectral(nEVP=48,coordinateKind="z");
            full = IMModes.geostrophicGeneralizedPotentialEnstrophyModes(background,1e-3,nModes=2,solver=solver);
            requested = IMModes.geostrophicGeneralizedPotentialEnstrophyModes(background,1e-3,nModes=0,solver=solver);
            selected = full.selectColumns([]);
            testCase.verifyEqual(requested.normalization,"generalizedPotentialEnstrophy");
            testCase.verifyEqual(requested.capabilities,selected.capabilities);
            testCase.verifyTrue(ismember(requested.normalization,requested.capabilities.normalizationNames));
            testCase.verifySize(requested.F([-1000;-500;0]),[3 0]);
            testCase.verifySize(requested.evaluate([-1000;-500;0],"Fz"),[3 0]);
            collection = IMBasisCollection({requested,selected});
            testCase.verifyEqual(collection.metadata(1).capabilities,collection.metadata(2).capabilities);
            testCase.verifySize(collection.evaluate([-1000;-500;0],variable="Fz"),[3 0 2]);
        end

        function knownIntrinsicSectorsDeclareTheirFirstFDerivative(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            barotropic = IMModes.hydrostaticFModes(background,nInternalModes=0);
            testCase.verifyEqual(barotropic.capabilities.derivativeOrders,{0:1,0:2});
            testCase.verifyEqual(barotropic.evaluate([-1000;0],"Fz"),zeros(2,1));
            energy = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf);
            endpoint = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1e-3,energyBoundary=energy);
            testCase.verifyEqual(endpoint.capabilities.derivativeOrders,{0:1,0:2});
            collection = IMBasisCollection({endpoint},basisIndex=[1 1]);
            testCase.verifySize(collection.evaluate([-1000;0],variable="Fz"),[2 0 2]);
            testCase.verifyError(@() collection.evaluate(-500,variable="Fz",derivativeOrder=1),"IMBasisCollection:UnsupportedDerivative");
        end
    end
end
