classdef IMModeBasisTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            originalPath = path;
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            testCase.addTeardown(@() path(originalPath));
        end
    end

    methods (Test)
        function numericalAndAnalyticalWavesShareTheDepthConvention(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            k = 1e-3;
            evp = IMInternalModes.waveModesAtWavenumber(N2=background.N2,zDomain=background.zDomain,k=k,rotationRate=background.rotationRate,latitude=background.latitude);
            numerical = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(evp,nModes=4);
            solution = background.analyticalSolution();
            analytical = solution.internalModes(evp,nModes=4);
            numericBasis = testCase.wrap(background,numerical,"waveEnergy",true,k);
            exactBasis = testCase.wrap(background,analytical,"waveEnergy",true,k);
            z = linspace(-1000,0,131).';
            for basis = [numericBasis,exactBasis]
                testCase.verifyEqual(basis.F(z),1000*basis.Gz(z),AbsTol=1e-13)
                testCase.verifyEqual(basis.evaluate(z,"F",1),1000*basis.Gzz(z),AbsTol=1e-12)
                for j = 1:4
                    value = integral(@(zz) waveEnergyDensity(basis,zz,j,k),-1000,0,AbsTol=1e-10)/1000;
                    testCase.verifyEqual(value,1,AbsTol=2e-10)
                end
            end
            testCase.verifyEqual(numericBasis.F(z),exactBasis.F(z),AbsTol=2e-9)
            testCase.verifyEqual(numericBasis.G(z),exactBasis.G(z),AbsTol=2e-9)
            testCase.verifyEqual(numericBasis.h,exactBasis.h,RelTol=2e-9)
        end

        function normalizationsScaleBothWaveProfilesTogether(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            evp = IMInternalModes.waveModesAtWavenumber(N2=background.N2,zDomain=background.zDomain,k=1e-3);
            solution = background.analyticalSolution();
            source = solution.internalModes(evp,nModes=3);
            rmsF = testCase.wrap(background,source,"depthRMSF",true,1e-3);
            rmsG = testCase.wrap(background,source,"depthRMSG",true,1e-3);
            z = linspace(-1000,0,201).';
            testCase.verifyEqual(trapz(z,rmsF.F(z).^2)/1000,ones(1,3),AbsTol=2e-13)
            testCase.verifyEqual(trapz(z,rmsG.G(z).^2)/1000,ones(1,3),AbsTol=2e-13)
            testCase.verifyEqual(rmsF.F(z).*rmsG.G(z),rmsG.F(z).*rmsF.G(z),AbsTol=1e-14)
            % Evaluating a copied source under another legacy normalization
            % cannot invalidate immutable facade normalization factors.
            expected = rmsF.F(z);
            source.normalization = "uMax";
            source.F(z);
            testCase.verifyEqual(rmsF.F(z),expected)
            testCase.verifyEqual(rmsF.normalization,"depthRMSF")
        end

        function analyticalSecondDerivativesRespectRotatingWaveEquation(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            k = 1e-3;
            evp = IMInternalModes.waveModesAtWavenumber(N2=background.N2,zDomain=background.zDomain,k=k,rotationRate=background.rotationRate,latitude=background.latitude);
            solution = background.analyticalSolution();
            source = solution.internalModes(evp,nModes=4);
            basis = testCase.wrap(background,source,"waveEnergy",true,k);
            z = linspace(-1000,0,131).';
            expected = (k^2-(background.N2(z)-background.f0^2)/background.g.*basis.eigenvalues).*basis.G(z);
            testCase.verifyEqual(basis.Gzz(z),expected,AbsTol=1e-16)
            testCase.verifyEqual(basis.evaluate(z,"Gz"),basis.Gz(z))
            testCase.verifyEqual(basis.evaluate(z,"Gzz"),basis.Gzz(z))
        end

        function signedMDAEnergyRetainsItsMetricSignature(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            g0 = -.04;
            gd = .03;
            evp = IMInternalModes.meanDensityAnomalyModes(N2=background.N2,zDomain=background.zDomain,g0=g0,gd=gd);
            source = IMSolverSpectral(nEVP=96,coordinateKind="z").solveEVP(evp,nModes=5);
            basis = testCase.wrap(background,source,"generalizedEnergy",false,0);
            [z,weights] = source.solver.nativeQuadratureRule(background.zDomain);
            G = basis.G(z);
            endpoints = basis.G(background.zDomain.');
            gram = G'*(weights.*(background.N2(z)/background.g).*G)+(g0/background.g)*(endpoints(2,:)'*endpoints(2,:))+(gd/background.g)*(endpoints(1,:)'*endpoints(1,:));
            testCase.verifyEqual(gram,diag(basis.normSignatures),AbsTol=2e-9)
            testCase.verifyTrue(any(basis.normSignatures < 0))
            testCase.verifyEqual(basis.F(0),zeros(1,5),AbsTol=1e-13)
            testCase.verifyEqual(basis.Gz(z),basis.evaluate(z,"G",1))
        end

        function generalizedEnstrophyUsesTheSameNamedPairing(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(N2=background.N2,zDomain=background.zDomain,k=1e-3,rotationRate=background.rotationRate,latitude=background.latitude);
            source = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(evp,nModes=4);
            basis = testCase.wrap(background,source,"generalizedPotentialEnstrophy",false,1e-3);
            z = linspace(-1000,0,31).';
            testCase.verifyEqual(basis.F(z),source.F(z,normalization="generalizedPotentialEnstrophy"),AbsTol=2e-13,RelTol=5e-14)
            testCase.verifyEqual(basis.G(z),source.G(z,normalization="generalizedPotentialEnstrophy"),AbsTol=2e-13,RelTol=5e-14)
            testCase.verifyEqual(basis.Gzz(z),-background.N2(z)/background.g.*basis.G(z)./basis.h,AbsTol=1e-13)
        end

        function neutralLabelsExposeOnlyTheirScientificChannels(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            source = IMNeutralLabelModesBasis([-1000 0],zDomain=background.zDomain,nModes=2);
            basis = testCase.wrap(background,source,"depthRMSG",false,0);
            z = [-1000;-700;-100;0];
            testCase.verifyEqual(basis.nColumns,3)
            testCase.verifyEqual(basis.G(z),source.G(z))
            testCase.verifyEqual(basis.Gz(z),source.Gz(z))
            testCase.verifyEqual(basis.Gzz(z),source.Gzz(z))
            testCase.verifyEqual(basis.APV(z),source.APV(z))
            testCase.verifyFalse(ismember("F",basis.capabilities))
            testCase.verifyError(@() basis.F(z),"IMModeBasis:UnavailableVariable")
            testCase.verifyError(@() basis.h,"IMModeBasis:UnavailableEquivalentDepth")
            testCase.verifyError(@() basis.eigenvalues,"IMModeBasis:UnavailableEigenvalues")
            testCase.verifyError(@() testCase.wrap(background,source,"depthRMSF",false,0),"IMModeBasis:UnavailableNormalization")
        end

        function scaledSelectedNeutralProjectionRecoversItsOwnCoefficients(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            source = IMNeutralLabelModesBasis([-800 -600],zDomain=background.zDomain,nModes=3);
            info = testCase.info(source);
            basis = IMModeBasis(source=source,background=background,family="neutralLabel",columns=[3 1],modeInfo=info([3 1],:),normalization="maxAbsG");
            z = linspace(-800,-600,4001).';
            weights = .05*ones(size(z));
            weights([1 end]) = .025;
            projection = basis.projectionOnGrid(z,weights);
            coefficients = [2.5;-.75];
            testCase.verifyEqual(projection.project(basis.APV(z)*coefficients),coefficients,AbsTol=2e-12)
            testCase.verifyEqual(projection.columnLabels,basis.columnLabels)
            originalG = source.G(z);
            testCase.verifyEqual(basis.G(z),originalG(:,[3 1])./sqrt(10),AbsTol=1e-14)
        end

        function constantNeutralLabelRequiresExternalAmplitudeInformation(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            source = IMNeutralLabelModesBasis([-1000 0],zDomain=background.zDomain,nModes=2);
            basis = testCase.wrap(background,source,"depthRMSG",false,0);
            z = linspace(-1000,0,1001).';
            weights = ones(size(z));
            weights([1 end]) = .5;
            testCase.verifyError(@() basis.projectionOnGrid(z,weights),"IMNeutralLabelModesBasis:UnobservableConstantLabel")
            projection = basis.projectionOnGrid(z,weights,columns=[2 3]);
            testCase.verifyEqual(projection.forwardMatrix*basis.APV(z),[zeros(2,1),eye(2)],AbsTol=2e-12)
        end

        function endpointResponsesKeepCanonicalScaling(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            k = 1e-4;
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=background.N2,zDomain=background.zDomain,k=k,rotationRate=background.rotationRate,latitude=background.latitude,surfaceBoundary=background.surfaceCondition);
            source = IMSolverSpectral(nEVP=64,coordinateKind="z").solveGeostrophicZeroAPVModes(problem);
            basis = testCase.wrap(background,source,"endpointResponse",false,k);
            z = linspace(-1000,0,31).';
            testCase.verifyEqual(basis.F(z),source.F(z))
            testCase.verifyEqual(basis.G(z),source.G(z))
            endpointF = basis.F(background.zDomain.');
            endpointG = basis.G(background.zDomain.');
            testCase.verifyEqual([endpointG(2,:)-endpointF(2,:);endpointG(1,:)],eye(2),AbsTol=2e-11)
            testCase.verifyEqual(basis.Gz(z),-background.g*k^2/background.f0^2*basis.F(z),AbsTol=1e-14)
        end

        function intrinsicAndEmptyBasesHaveTheSameShapeContract(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            info = table("critical:1","critical",1,VariableNames=["columnLabel","branch","branchIndex"]);
            evaluator = @(variable,z,order) intrinsicWave(variable,z,order);
            basis = IMModeBasis(background=background,family="zeroDepthWave",modeInfo=info,normalization="waveEnergy",waveConvention=true,k=1e-3,h=0,intrinsicEvaluator=evaluator,variables=["F","G"]);
            z = linspace(-1000,0,31).';
            testCase.verifyEqual(basis.F(z),1000*basis.Gz(z),AbsTol=1e-14)
            testCase.verifyEqual(basis.h,0)
            emptyBasis = IMModeBasis(background=background,family="wave",normalization="waveEnergy",waveConvention=true,variables=["F","G"]);
            testCase.verifySize(emptyBasis.G(z),[31 0])
            testCase.verifySize(emptyBasis.F(z),[31 0])
            testCase.verifyEqual(emptyBasis.nColumns,0)
            testCase.verifyError(@() IMModeBasis(background=background,family="null",modeInfo=info,normalization="depthRMSF",intrinsicEvaluator=@(~,zz,~) zeros(size(zz))),"IMModeBasis:NullNormalization")
        end

        function maximumNormalizationResolvesInteriorExtrema(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            evp = IMInternalModes.waveModesAtWavenumber(N2=background.N2,zDomain=background.zDomain,k=1e-3);
            numerical = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(evp,nModes=3);
            solution = background.analyticalSolution();
            analytical = solution.internalModes(evp,nModes=3);
            for source = {numerical,analytical}
                maxF = testCase.wrap(background,source{1},"maxAbsF",true,1e-3);
                maxG = testCase.wrap(background,source{1},"maxAbsG",true,1e-3);
                testCase.verifyEqual(abs(maxF.F(0)),ones(1,3),AbsTol=2e-12)
                for j = 1:3
                    values = maxG.G(-1000+1000/(2*j));
                    testCase.verifyEqual(abs(values(j)),1,AbsTol=2e-12)
                end
            end
        end

        function intrinsicNormalizationSplitsExactSupportIntervals(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            intervals = [-812 -631;-111 -34];
            info = table(["critical:1";"critical:2"],repmat("critical",2,1),ones(2,1),VariableNames=["columnLabel","branch","branchIndex"]);
            evaluator = @(variable,z,order) compactProfiles(variable,z,order,intervals);
            basis = IMModeBasis(background=background,family="zeroDepthWave",modeInfo=info,normalization="depthRMSG",waveConvention=true,h=[0 0],intrinsicEvaluator=evaluator,integrationIntervals=intervals);
            peaks = mean(intervals,2);
            expected = diag(sqrt(2000./diff(intervals,1,2)));
            testCase.verifyEqual(basis.G(peaks),expected,AbsTol=2e-13)
            testCase.verifyEqual(basis.F(peaks),1000*basis.Gz(peaks),AbsTol=1e-13)
        end
    end

    methods (Static, Access = private)
        function basis = wrap(background,source,normalization,waveConvention,k)
            basis = IMModeBasis(source=source,background=background,family="test",modeInfo=IMModeBasisTests.info(source),normalization=normalization,waveConvention=waveConvention,k=k);
        end

        function info = info(source)
            if isprop(source,"columnLabels")
                labels = source.columnLabels;
            elseif isprop(source,"endpoints")
                labels = source.endpoints;
            else
                labels = "mode:"+string(source.modeNumber);
            end
            n = numel(labels);
            info = table(labels(:),repmat("test",n,1),(1:n).',VariableNames=["columnLabel","branch","branchIndex"]);
        end
    end
end

function values = waveEnergyDensity(basis,z,j,k)
F = basis.F(z(:));
G = basis.G(z(:));
values = reshape(F(:,j).^2+(k*diff(basis.zDomain))^2*G(:,j).^2,size(z));
end

function values = intrinsicWave(variable,z,order)
m = pi/1000;
phase = (z+1000)*m;
if variable == "F"
    order = order+1;
    factor = 1000;
else
    factor = 1;
end
values = factor*m^order*sin(phase+order*pi/2);
end

function values = compactProfiles(variable,z,order,intervals)
if variable == "F"
    order = order+1;
    factor = 1000;
else
    factor = 1;
end
values = zeros(numel(z),size(intervals,1));
for j = 1:size(intervals,1)
    inside = z >= intervals(j,1) & z <= intervals(j,2);
    m = pi/diff(intervals(j,:));
    values(inside,j) = factor*m^order*sin((z(inside)-intervals(j,1))*m+order*pi/2);
end
end
