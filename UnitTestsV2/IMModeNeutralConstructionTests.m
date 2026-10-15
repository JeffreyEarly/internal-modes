classdef IMModeNeutralConstructionTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function fullyNeutralRigidHydrostaticFKeepsItsBarotropicMode(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = IMModes.hydrostaticFModes(background,solver=IMSolverSpectral(nEVP=32),nInternalModes=0);
            testCase.verifyEqual(basis.nColumns,1);
            testCase.verifyEqual(basis.modeInfo.branch,"barotropic");
            testCase.verifyEqual(basis.eigenvalues,0);
            testCase.verifyEqual(basis.h,Inf);
            testCase.verifyEqual(basis.normalization,"depthRMSF");
            z = linspace(-1000,0,13).';
            testCase.verifyEqual(basis.F(z),ones(size(z)),AbsTol=2e-13);
            testCase.verifyEqual(basis.G(z),zeros(size(z)),AbsTol=2e-13);
        end

        function rigidHydrostaticFUsesTheSameRegularAnalyticalAndNumericalPair(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            exact = IMModes.hydrostaticFModes(background,solver=IMSolverAnalytical(),nInternalModes=3);
            numerical = IMModes.hydrostaticFModes(background,solver=IMSolverSpectral(nEVP=64),nInternalModes=3);
            testCase.verifyEqual(numerical.nColumns,4);
            testCase.verifyEqual(numerical.modeInfo.branch,["barotropic";"positive";"positive";"positive"]);
            testCase.verifyEqual(numerical.columnLabels,exact.columnLabels);
            testCase.verifyEqual(numerical.h(2:end),exact.h(2:end),RelTol=2e-8);
            z = linspace(-1000,0,31).';
            testCase.verifyEqual(numerical.F(z),exact.F(z),AbsTol=2e-9);
            testCase.verifyEqual(numerical.G(z),exact.G(z),AbsTol=2e-9);
        end

        function balancedZeroCountsRetainAdmissibleSpecialSectors(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            inactive = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48)}
                apv = IMModes.geostrophicAPVModes(background,solver=solver{1},energyBoundary=inactive,nPositiveModes=0);
                hydrostatic = IMModes.hydrostaticFModes(background,solver=solver{1},nInternalModes=0);
                testCase.verifyEqual(apv.nColumns,1);
                testCase.verifyEqual(hydrostatic.nColumns,1);
                testCase.verifyEqual(apv.eigenvalues,0);
                testCase.verifyEqual(hydrostatic.eigenvalues,0);
                testCase.verifyEqual(apv.modeInfo.branch,"zeroEigenvalue");
                testCase.verifyEqual(hydrostatic.modeInfo.branch,"barotropic");
            end
            neumann = IMGeneralizedEnergyBoundary(g0=0,gd=0);
            mda = IMModes.meanDensityAnomalyModes(background,solver=IMSolverSpectral(nEVP=48),energyBoundary=neumann,nPositiveModes=0);
            testCase.verifyEqual(mda.nColumns,1);
            testCase.verifyEqual(mda.modeInfo.branch,"constant");
            testCase.verifyEqual(mda.eigenvalues,0);
            G = mda.G([-1000;-500;0]);
            testCase.verifyEqual(G,repmat(G(1,:),3,1),AbsTol=2e-12);
        end

        function zeroCountsKeepEmptyEVPSpectraAvailable(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48)}
                wave = IMModes.waveModesAtWavenumber(background,2e-4,solver=solver{1},nInternalModes=0,nSubinertialModes=0);
                inertial = IMModes.inertialModes(background,solver=solver{1},nPositiveModes=0,nNegativeModes=0);
                hydrostatic = IMModes.hydrostaticGModes(background,solver=solver{1},nInternalModes=0);
                for basis = {wave,inertial,hydrostatic}
                    testCase.verifyEqual(basis{1}.nColumns,0);
                    testCase.verifyEmpty(basis{1}.h);
                    testCase.verifyEmpty(basis{1}.eigenvalues);
                    testCase.verifySize(basis{1}.F([-1000;0]),[2 0]);
                    testCase.verifySize(basis{1}.G([-1000;0]),[2 0]);
                end
            end
        end

        function inactiveZeroAPVEndpointsHaveNoFictitiousSpectrum(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            inactive = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=32)}
                basis = IMModes.geostrophicZeroAPVModesAtWavenumber(background,2e-4,solver=solver{1},energyBoundary=inactive);
                testCase.verifyEqual(basis.nColumns,0);
                testCase.verifyEmpty(basis.columnLabels);
                testCase.verifySize(basis.F([-1000;0]),[2 0]);
                testCase.verifySize(basis.G([-1000;0]),[2 0]);
                testCase.verifyError(@() basis.h,"IMModeBasis:UnavailableEquivalentDepth");
                testCase.verifyError(@() basis.eigenvalues,"IMModeBasis:UnavailableEigenvalues");
            end
        end

        function numericalAndExactEndpointsKeepTheSameAutomaticResponses(testCase)
            background = IMModeBackground.constantStratification(N0=.005,zDomain=[-1000 0]);
            z = linspace(-1000,0,41).';
            for energy = {IMGeneralizedEnergyBoundary(g0=0,gd=0),IMGeneralizedEnergyBoundary(g0=-.02,gd=Inf),IMGeneralizedEnergyBoundary(g0=Inf,gd=.02)}
                exact = IMModes.geostrophicZeroAPVModesAtWavenumber(background,2e-4,solver=IMSolverAnalytical(),energyBoundary=energy{1});
                numerical = IMModes.geostrophicZeroAPVModesAtWavenumber(background,2e-4,solver=IMSolverSpectral(nEVP=64),energyBoundary=energy{1});
                testCase.verifyEqual(numerical.columnLabels,exact.columnLabels);
                testCase.verifyEqual(numerical.nColumns,numel(energy{1}.activeEndpoints));
                testCase.verifyEqual(numerical.F(z),exact.F(z),AbsTol=2e-9);
                testCase.verifyEqual(numerical.G(z),exact.G(z),AbsTol=2e-9);
            end
        end
    end
end
