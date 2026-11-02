classdef IMHorizontalIdentityTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function fixedFrequencyRootsSatisfyRigidLidDispersion(testCase)
            % The exact vertical wavenumbers give a reference independent
            % of the returned equivalent-depth and eigenvalue arrays.
            N0 = .01; D = 1000;
            for latitude = [-33 33]
                background = IMModeBackground.constantStratification(N0=N0,zDomain=[-D 0],latitude=latitude,surfaceCondition="rigidLid");
                for omega = [.5*abs(background.f0) .005]
                    for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48,coordinateKind="z")}
                        if omega > abs(background.f0)
                            basis = IMModes.waveModesAtFrequency(background,omega,solver=solver{1},nPropagatingModes=3,nEvanescentModes=0);
                        else
                            basis = IMModes.waveModesAtFrequency(background,omega,solver=solver{1},nPropagatingModes=0,nEvanescentModes=3);
                        end
                        expected = (omega^2-background.f0^2)*(pi*(1:3)/D).^2/(N0^2-omega^2);
                        testCase.verifyEqual(basis.horizontalKind,"perColumn");
                        testCase.verifyEmpty(basis.k);
                        testCase.verifyEqual(basis.kSquared,expected,RelTol=3e-11);
                        testCase.verifyGreaterThan(basis.h,zeros(1,3));
                        if omega < abs(background.f0)
                            % Positive h does not imply a propagating root.
                            testCase.verifyLessThan(basis.kSquared,zeros(1,3));
                        end
                        testCase.verifyEqual(basis.frequency,repmat(omega,1,3));
                        testCase.verifyEqual(basis.scientificSpecification.resolved.horizontal,struct(kind="perColumn",k=zeros(1,0),kSquared=basis.kSquared));
                    end
                end
            end
        end

        function selectingRootsPreservesDispersionAndScientificIdentity(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = IMModes.waveModesAtFrequency(background,.005,solver=IMSolverAnalytical(),nPropagatingModes=4);
            selected = basis.selectColumns([4 1 3]).selectColumns([3 1]);
            columns = [3 4];
            testCase.verifyEqual(selected.kSquared,basis.kSquared(columns));
            testCase.verifyEqual(selected.h,basis.h(columns));
            testCase.verifyEqual(selected.eigenvalues,basis.eigenvalues(columns));
            testCase.verifyEqual(selected.columnLabels,basis.columnLabels(columns));
            testCase.verifyEqual(selected.constructionRecord.selection.initialColumnIndices,columns);
            testCase.verifyEqual(selected.constructionRecord.requested,basis.constructionRecord.requested);
            z = [-900;-510;-50];
            profiles = basis.F(z);
            testCase.verifyEqual(selected.F(z),profiles(:,columns));
            empty = selected.selectColumns([]);
            testCase.verifyEqual(empty.horizontalKind,"perColumn");
            testCase.verifyEmpty(empty.kSquared);
            testCase.verifyEmpty(empty.scientificSpecification.resolved.horizontal.kSquared);
            testCase.verifyTrue(ismember("kSquared",empty.capabilities.spectrumNames));
            testCase.verifySize(empty.F(z),[3 0]);
        end

        function independentFamiliesAcceptExternallyAssignedPages(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            energy = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf);
            neutralBackground = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            bases = {IMModes.geostrophicAPVModes(background,solver=IMSolverAnalytical(),energyBoundary=energy,nPositiveModes=2), ...
                IMModes.meanDensityAnomalyModes(background,solver=IMSolverSpectral(nEVP=48),energyBoundary=energy,nPositiveModes=2), ...
                IMModes.hydrostaticFModes(background,solver=IMSolverAnalytical(),nInternalModes=2), ...
                IMModes.hydrostaticGModes(background,solver=IMSolverAnalytical(),nInternalModes=2), ...
                IMModes.neutralLabelModes(neutralBackground,neutralIntervals=[-1000 0],nModesPerInterval=2)};
            for basisCell = bases
                basis = basisCell{1};
                testCase.verifyEqual(basis.horizontalKind,"independent");
                testCase.verifyEmpty(basis.k);
                testCase.verifyEmpty(basis.kSquared);
                collection = IMBasisCollection({basis},kappa=[1e-3 3e-3 1e-3],basisIndex=[1 1 1]);
                testCase.verifyEqual(collection.kappa,[1e-3 3e-3 1e-3]);
                testCase.verifyEqual(collection.metadata.horizontalIdentity.kind,"independent");
                testCase.verifyEmpty(collection.metadata.horizontalIdentity.k);
                testCase.verifyEqual(collection.metadata.scientificSpecification,basis.scientificSpecification);
                testCase.verifyEmpty(IMBasisCollection({basis}).kappa);
            end
        end

        function solvedWavenumbersCannotBeReassigned(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            wave = IMModes.waveModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),nInternalModes=2);
            inertial = IMModes.inertialModes(background,solver=IMSolverAnalytical(),nPositiveModes=2,nNegativeModes=0);
            endpoint = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),energyBoundary=IMGeneralizedEnergyBoundary(g0=0,gd=0));
            enstrophy = IMModes.geostrophicGeneralizedPotentialEnstrophyModes(background,1e-3,solver=IMSolverSpectral(nEVP=48),nModes=2);
            for basisCell = {wave,inertial,endpoint,enstrophy}
                basis = basisCell{1};
                testCase.verifyEqual(basis.horizontalKind,"fixedWavenumber");
                testCase.verifyEmpty(basis.kSquared);
                testCase.verifyEqual(IMBasisCollection({basis}).kappa,basis.k);
                testCase.verifyError(@() IMBasisCollection({basis},kappa=basis.k+1e-3),"IMBasisCollection:KappaMismatch");
            end
            testCase.verifyEqual(inertial.k,0);
            testCase.verifyEqual(IMBasisCollection({inertial},kappa=0).kappa,0);
        end

        function perColumnRootsRemainDistinctFromCollectionPages(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = IMModes.waveModesAtFrequency(background,.005,solver=IMSolverAnalytical(),nPropagatingModes=3);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            testCase.verifyEmpty(collection.kappa);
            testCase.verifyEqual(collection.metadata.horizontalIdentity.kSquared,basis.kSquared);
            testCase.verifyError(@() IMBasisCollection({basis},kappa=1e-3),"IMBasisCollection:UnsupportedHorizontalPage");
            empty = IMModes.waveModesAtFrequency(background,.005,solver=IMSolverAnalytical(),nPropagatingModes=0,nEvanescentModes=0);
            testCase.verifyEqual(empty.horizontalKind,"perColumn");
            testCase.verifyEmpty(empty.k);
            testCase.verifyEmpty(empty.kSquared);
            testCase.verifyEqual(empty.capabilities.spectrumNames,["h","eigenvalues","kSquared","frequency"]);
            testCase.verifyEqual(empty.capabilities,basis.selectColumns([]).capabilities);
            testCase.verifyError(@() IMBasisCollection({empty},kappa=1e-3),"IMBasisCollection:UnsupportedHorizontalPage");
        end

        function directWrappersMustDeclareTheirHorizontalMeaning(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            info = table("profile:1","intrinsic",1,VariableNames=["columnLabel","branch","branchIndex"]);
            basis = IMModeBasis(background=background,family="intrinsic",modeInfo=info,normalization="depthRMSF",intrinsicEvaluator=@(~,z,~) ones(size(z)),capabilities="F");
            testCase.verifyEqual(basis.horizontalKind,"unspecified");
            testCase.verifyEmpty(basis.k);
            testCase.verifyEmpty(IMBasisCollection({basis}).kappa);
            testCase.verifyError(@() IMBasisCollection({basis},kappa=1e-3),"IMBasisCollection:UnsupportedHorizontalPage");
            testCase.verifyError(@() IMModeBasis(background=background,family="intrinsic",horizontalKind="independent",k=1e-3),"IMModeBasis:InvalidHorizontalIdentity");
            testCase.verifyError(@() IMModeBasis(background=background,family="intrinsic",horizontalKind="fixedWavenumber"),"IMModeBasis:InvalidHorizontalIdentity");
            testCase.verifyError(@() IMModeBasis(background=background,family="intrinsic",horizontalKind="perColumn",kSquared=[1 -1]),"IMModeBasis:InvalidHorizontalIdentity");
        end
    end
end
