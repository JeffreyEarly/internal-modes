classdef IMModeCreationTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function bulkRequestsPreservePagesAndIndependentCounts(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            collection = IMModes.waveModesAtWavenumbers(background,[2e-4 1e-4 2e-4],solver=IMSolverAnalytical(),nInternalModes=[2 3 2]);
            testCase.verifyEqual(collection.kappa,[2e-4 1e-4 2e-4])
            testCase.verifyEqual(collection.basisIndex,[1 2 1])
            testCase.verifyEqual(numel(collection.bases),2)
            testCase.verifyEqual(collection.bases{1}.nColumns,3)
            testCase.verifyEqual(collection.bases{2}.nColumns,4)
            z = linspace(-1000,0,21).';
            F = collection.evaluate(z,variable="F",pages=[3 1]);
            testCase.verifyEqual(F(:,:,1),F(:,:,2))
            testCase.verifyError(@() IMModes.waveModesAtWavenumbers(background,[1e-4 1e-4],nInternalModes=[2 3]),"IMModes:ConflictingModeCounts")
        end

        function intervalDeclarationsAndNormalizationContractsAreValidated(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            energy = IMGeneralizedEnergyBoundary(g0=0,gd=0);
            testCase.verifyError(@() IMModes.neutralLabelModes(background,neutralIntervals=[-200 0],energyBoundary=energy),"IMModes:InconsistentInterval")
            testCase.verifyError(@() IMModes.zeroDepthWaveModesAtWavenumber(background,1e-4,criticalIntervals=[-200 0]),"IMModes:InconsistentInterval")
            testCase.verifyError(@() IMModes.waveModesAtFrequency(background,.002,normalization="waveEnergy"),"IMModes:UnavailableNormalization")
            testCase.verifyError(@() IMModes.waveModesAtFrequency(background,abs(background.f0)),"IMModes:InertialFrequency")
            testCase.verifyError(@() IMModes.meanDensityAnomalyModes(background,solver=IMSolverAnalytical(),energyBoundary=energy),"IMModes:UnsupportedAnalyticalFamily")
        end

        function criticalRigidColumnHasAnEmptyRegularWaveInventory(testCase)
            f0 = 2*7.2921e-5*sind(33);
            background = IMModeBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = IMModes.waveModesAtWavenumber(background,1e-4,nInternalModes=0,nSubinertialModes=0);
            testCase.verifyEqual(basis.nColumns,0)
            testCase.verifyEqual(basis.h,zeros(1,0))
            testCase.verifyError(@() IMModes.waveModesAtWavenumber(background,1e-4,nInternalModes=1),"IMModes:InsufficientModes")
        end

        function wavesKeepSeparateLaddersAndAutomaticSurfaceMode(testCase)
            N2 = @(z) 1e-4*max(min((-z-100)/100,1),0).^2;
            background = IMModeBackground(N2=N2,zDomain=[-1000 0]);
            k = 2e-4;
            basis = IMModes.waveModesAtWavenumber(background,k,solver=IMSolverSpectral(nEVP=160,coordinateKind="z"),nInternalModes=3,nSubinertialModes=2);
            testCase.verifyClass(basis,"IMModeBasis")
            testCase.verifyEqual(basis.nColumns,6)
            testCase.verifyEqual(basis.modeInfo.branch,["surface";"internal";"internal";"internal";"subinertial";"subinertial"])
            testCase.verifyEqual(basis.modeInfo.branchIndex,[1;1;2;3;1;2])
            testCase.verifyEqual(string(basis.modeInfo.columnLabel),basis.columnLabels(:))
            testCase.verifyEqual(basis.normalization,"waveEnergy")
            testCase.verifyGreaterThan(basis.h(1:4),zeros(1,4))
            testCase.verifyLessThan(basis.h(5:6),zeros(1,2))
            testCase.verifyLessThan(diff(basis.frequency(2:4)),zeros(1,2))
            testCase.verifyGreaterThan(diff(basis.frequency(5:6)),0)
            testCase.verifyGreaterThan(basis.frequency(2:4),abs(background.f0)*ones(1,3))
            testCase.verifyLessThan(basis.frequency(5:6),abs(background.f0)*ones(1,2))
            z = linspace(-1000,0,137).';
            testCase.verifyEqual(basis.F(z),background.D*basis.Gz(z),RelTol=2e-11,AbsTol=2e-11)
            for j = 1:basis.nColumns
                norm = integral(@(z) waveNormIntegrand(basis,z,j,k,background.D),-1000,0,AbsTol=1e-8,RelTol=1e-10)/background.D;
                testCase.verifyEqual(norm,1,AbsTol=2e-8)
            end
        end

        function zeroWaveCountsReturnOnlyAdmissibleSpecialModes(testCase)
            free = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            rigid = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            surface = IMModes.waveModesAtWavenumber(free,2e-4,solver=IMSolverAnalytical(),nInternalModes=0,nSubinertialModes=0);
            empty = IMModes.waveModesAtWavenumber(rigid,2e-4,solver=IMSolverAnalytical(),nInternalModes=0,nSubinertialModes=0);
            testCase.verifyEqual(surface.nColumns,1)
            testCase.verifyEqual(surface.modeInfo.branch,"surface")
            testCase.verifyEqual(empty.nColumns,0)
            testCase.verifySize(empty.G([-1000;0]),[2 0])
        end

        function analyticalAndNumericalWavesShareIdentityAndConvention(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],latitude=-30);
            k = 2e-4;
            exact = IMModes.waveModesAtWavenumber(background,k,solver=IMSolverAnalytical(),nInternalModes=4,nSubinertialModes=0);
            numerical = IMModes.waveModesAtWavenumber(background,k,solver=IMSolverSpectral(nEVP=64,coordinateKind="z"),nInternalModes=4,nSubinertialModes=0);
            testCase.verifyEqual(numerical.columnLabels,exact.columnLabels)
            testCase.verifyEqual(numerical.modeInfo.branch,exact.modeInfo.branch)
            testCase.verifyEqual(numerical.modeInfo.branchIndex,exact.modeInfo.branchIndex)
            testCase.verifyEqual(numerical.h,exact.h,RelTol=2e-7)
            z = linspace(-1000,0,101).';
            testCase.verifyEqual(exact.F(z),1000*exact.Gz(z),AbsTol=2e-12)
            testCase.verifyEqual(numerical.F(z),exact.F(z),RelTol=2e-6,AbsTol=2e-7)
            testCase.verifyEqual(numerical.G(z),exact.G(z),RelTol=2e-6,AbsTol=2e-7)
        end

        function waveNormalizationsShareNamesAndScalePairedProfiles(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            k = 2e-4;
            z = linspace(-1000,0,20001).';
            names = ["depthRMSF","depthRMSG","maxAbsF","maxAbsG"];
            for name = names
                basis = IMModes.waveModesAtWavenumber(background,k,solver=IMSolverAnalytical(),nInternalModes=2,nSubinertialModes=0,normalization=name);
                testCase.verifyEqual(basis.normalization,name)
                F = basis.F(z);
                G = basis.G(z);
                testCase.verifyEqual(F,1000*basis.Gz(z),RelTol=2e-12,AbsTol=2e-12)
                switch name
                    case "depthRMSF"
                        norms = trapz(z,abs(F).^2)/1000;
                    case "depthRMSG"
                        norms = trapz(z,abs(G).^2)/1000;
                    case "maxAbsF"
                        norms = max(abs(F),[],1);
                    case "maxAbsG"
                        norms = max(abs(G),[],1);
                end
                testCase.verifyEqual(norms,ones(1,basis.nColumns),AbsTol=2e-7)
            end
        end

        function apvAlwaysIncludesFiniteNegativeAndZeroSectors(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            negativeBoundary = IMGeneralizedEnergyBoundary(g0=-.04,gd=0);
            zeroBoundary = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=64,coordinateKind="z")}
                negative = IMModes.geostrophicAPVModes(background,solver=solver{1},energyBoundary=negativeBoundary,nPositiveModes=0);
                zero = IMModes.geostrophicAPVModes(background,solver=solver{1},energyBoundary=zeroBoundary,nPositiveModes=0);
                testCase.verifyEqual(negative.nColumns,1)
                testCase.verifyLessThan(negative.eigenvalues,0)
                testCase.verifyEqual(negative.modeInfo.branch,"negative")
                testCase.verifyEqual(zero.nColumns,1)
                testCase.verifyEqual(zero.eigenvalues,0)
                testCase.verifyEqual(zero.h,Inf)
                testCase.verifyEqual(zero.modeInfo.branch,"zeroEigenvalue")
                positive = IMModes.geostrophicAPVModes(background,solver=solver{1},energyBoundary=negativeBoundary,nPositiveModes=3);
                testCase.verifyEqual(positive.nColumns,4)
                testCase.verifyEqual(nnz(positive.eigenvalues > 0),3)
                testCase.verifyEqual(positive.normalization,"depthRMSF")
                z = linspace(-1000,0,10001).';
                testCase.verifyEqual(trapz(z,abs(positive.F(z)).^2)/1000,ones(1,4),AbsTol=1e-7)
            end
        end

        function meanDensityConstantIsAutomatic(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            energyBoundary = IMGeneralizedEnergyBoundary(g0=.02,gd=-.01);
            solver = IMSolverSpectral(nEVP=64,coordinateKind="z");
            basis = IMModes.meanDensityAnomalyModes(background,solver=solver,energyBoundary=energyBoundary,nPositiveModes=0);
            testCase.verifyGreaterThanOrEqual(basis.nColumns,1)
            zero = basis.eigenvalues == 0;
            testCase.verifyEqual(nnz(zero),1)
            testCase.verifyEqual(nnz(basis.eigenvalues > 0),0)
            testCase.verifyEqual(basis.normalization,"generalizedEnergy")
            z = linspace(-1000,0,101).';
            G = basis.G(z);
            testCase.verifyEqual(G(:,zero),repmat(G(1,zero),numel(z),1),AbsTol=2e-10)
            regular = IMModes.meanDensityAnomalyModes(background,solver=solver,energyBoundary=energyBoundary,nPositiveModes=3);
            testCase.verifyEqual(nnz(regular.eigenvalues > 0),3)
            testCase.verifyEqual(nnz(regular.eigenvalues <= 0),basis.nColumns)
        end

        function zeroAPVReturnsEveryFiniteEndpoint(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            solver = IMSolverSpectral(nEVP=64,coordinateKind="z");
            both = IMModes.geostrophicZeroAPVModesAtWavenumber(background,2e-4,solver=solver,energyBoundary=IMGeneralizedEnergyBoundary(g0=0,gd=-.02));
            surface = IMModes.geostrophicZeroAPVModesAtWavenumber(background,2e-4,solver=solver,energyBoundary=IMGeneralizedEnergyBoundary(g0=.02,gd=Inf));
            bottom = IMModes.geostrophicZeroAPVModesAtWavenumber(background,2e-4,solver=solver,energyBoundary=IMGeneralizedEnergyBoundary(g0=Inf,gd=0));
            testCase.verifyEqual(both.nColumns,2)
            testCase.verifyEqual(surface.nColumns,1)
            testCase.verifyEqual(bottom.nColumns,1)
            testCase.verifyEqual(both.normalization,"endpointResponse")
            testCase.verifyEqual(both.columnLabels,[surface.columnLabels bottom.columnLabels])
            z = linspace(-1000,0,37).';
            testCase.verifyEqual(both.F(z),[surface.F(z) bottom.F(z)],RelTol=2e-10,AbsTol=2e-10)
        end

        function neutralLabelsKeepTheirConstantAndUnobservableAmplitude(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            energyBoundary = IMGeneralizedEnergyBoundary(g0=0,gd=0);
            basis = IMModes.neutralLabelModes(background,neutralIntervals=[-1000 0],energyBoundary=energyBoundary,nModesPerInterval=0);
            testCase.verifyEqual(basis.nColumns,1)
            testCase.verifyEqual(basis.normalization,"depthRMSG")
            z = linspace(-1000,0,11).';
            testCase.verifyEqual(basis.G(z),ones(size(z)))
            testCase.verifyEqual(basis.Gz(z),zeros(size(z)))
            testCase.verifyEqual(basis.APV(z),zeros(size(z)))
            testCase.verifyError(@() basis.projectionOnGrid(z,ones(size(z))),"IMNeutralLabelModesBasis:UnobservableConstantLabel")
            regular = IMModes.neutralLabelModes(background,neutralIntervals=[-1000 0],energyBoundary=energyBoundary,nModesPerInterval=3);
            testCase.verifyEqual(regular.nColumns,4)
            testCase.verifyEqual(regular.columnLabels(1),basis.columnLabels)
        end

        function hydrostaticFactoriesKeepExternalModeOutsideInternalCount(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=64,coordinateKind="z")}
                F = IMModes.hydrostaticFModes(background,solver=solver{1},nInternalModes=2);
                G = IMModes.hydrostaticGModes(background,solver=solver{1},nInternalModes=2);
                testCase.verifyEqual(F.nColumns,3)
                testCase.verifyEqual(G.nColumns,3)
                testCase.verifyEqual(F.h,G.h,RelTol=2e-7)
                testCase.verifyEqual(IMModes.hydrostaticGModes(background,solver=solver{1},nInternalModes=0).nColumns,1)
            end
        end

        function inertialModesExposeAuxiliaryLadders(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            basis = IMModes.inertialModes(background,solver=IMSolverAnalytical(),nPositiveModes=3,nNegativeModes=0);
            testCase.verifyEqual(basis.nColumns,4)
            testCase.verifyEqual(basis.normalization,"depthRMSF")
            testCase.verifyEqual(basis.frequency,abs(background.f0)*ones(1,4),AbsTol=2e-15)
            z = linspace(-1000,0,101).';
            testCase.verifyEqual(basis.F(z),1000*basis.Gz(z),AbsTol=2e-12)
            testCase.verifyEqual(basis.evaluate(z,"F",1),1000*basis.Gzz(z),AbsTol=2e-12)
        end

        function generalizedEnstrophyFactoryReturnsRequestedLadder(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            basis = IMModes.geostrophicGeneralizedPotentialEnstrophyModes(background,2e-4,solver=IMSolverSpectral(nEVP=64,coordinateKind="z"),nModes=3);
            testCase.verifyEqual(basis.nColumns,3)
            testCase.verifyEqual(basis.normalization,"generalizedPotentialEnstrophy")
            testCase.verifyGreaterThan(basis.eigenvalues,zeros(1,3))
        end

        function fixedFrequencyWavesUseTheSameProfileConvention(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            omega = .002;
            basis = IMModes.waveModesAtFrequency(background,omega,solver=IMSolverAnalytical(),nPropagatingModes=3,nEvanescentModes=0);
            testCase.verifyEqual(basis.nColumns,4)
            testCase.verifyEqual(basis.normalization,"depthRMSF")
            testCase.verifyEqual(basis.frequency,omega*ones(1,4),AbsTol=2e-15)
            z = linspace(-1000,0,10001).';
            F = basis.F(z);
            testCase.verifyEqual(F,1000*basis.Gz(z),AbsTol=2e-12)
            testCase.verifyEqual(trapz(z,abs(F).^2)/1000,ones(1,4),AbsTol=1e-7)
        end

        function exactZeroDepthWavesHaveFiniteRequestedAuxiliaryBasis(testCase)
            f0 = 2*7.2921e-5*sind(33);
            background = IMModeBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0]);
            basis = IMModes.zeroDepthWaveModesAtWavenumber(background,2e-4,criticalIntervals=[-1000 0],nModesPerInterval=2);
            testCase.verifyEqual(basis.nColumns,2)
            testCase.verifyEqual(basis.h,zeros(1,2))
            testCase.verifyEqual(basis.frequency,abs(f0)*ones(1,2),AbsTol=2e-15)
            z = linspace(-1000,0,101).';
            testCase.verifyEqual(basis.F(z),1000*basis.Gz(z),AbsTol=2e-12)
            testCase.verifyEqual(basis.evaluate(z,"F",1),1000*basis.Gzz(z),AbsTol=2e-12)
            empty = IMModes.zeroDepthWaveModesAtWavenumber(background,2e-4,criticalIntervals=[-1000 0],nModesPerInterval=0);
            testCase.verifyEqual(empty.nColumns,0)
            testCase.verifyEqual(empty.h,zeros(1,0))
        end

        function numericalRequestsAreExactAndRejectUnavailableSector(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            small = IMSolverSpectral(nEVP=16,coordinateKind="z");
            testCase.verifyError(@() IMModes.waveModesAtWavenumber(background,2e-4,solver=small,nInternalModes=30,nSubinertialModes=0),"IMModes:InsufficientModes")
            testCase.verifyError(@() IMModes.waveModesAtWavenumber(background,2e-4,solver=IMSolverAnalytical(),nInternalModes=0,nSubinertialModes=1),"IMModes:InsufficientModes")
            testCase.verifyError(@() IMModes.waveModesAtWavenumber(background,2e-4,nMixedLayerModes=1),"MATLAB:TooManyInputs")
            testCase.verifyError(@() IMModes.geostrophicAPVModes(background,nPositiveModes=2,negativeModes="exclude"),"MATLAB:TooManyInputs")
        end
    end
end

function values = waveNormIntegrand(basis,z,column,k,D)
F = basis.F(z(:));
G = basis.G(z(:));
values = reshape(abs(F(:,column)).^2+(k*D)^2*abs(G(:,column)).^2,size(z));
end
