classdef IMEigenvalueProblemRefactorTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.originalPath = path;
            addpath(repoRoot);
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function spectralWaveAtWavenumberAssemblyMatchesPhysicalStrongForm(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            k = 1e-4;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=f0, g=g);

            [A, B] = evp.assemble(solver);
            expectedA = solver.Txx - k*k*solver.T;
            expectedB = diag((f0*f0 - N2(solver.zNative))/g)*solver.T;

            interiorRows = 2:(nEVP-1);
            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(A(1,:), solver.T(1,:), AbsTol=1e-12)
            testCase.verifyEqual(B(1,:), zeros(1,nEVP), AbsTol=1e-12)
            testCase.verifyEqual(A(end,:), solver.T(end,:), AbsTol=1e-12)
            testCase.verifyEqual(B(end,:), zeros(1,nEVP), AbsTol=1e-12)
        end

        function spectralWaveAtFrequencyAssemblyMatchesPhysicalStrongForm(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            omega = 0.8*5.2e-3;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = IMEigenvalueProblem.waveModesAtFrequency(omega=omega, f0=f0, g=g);

            [A, B] = evp.assemble(solver);
            expectedA = solver.Txx;
            expectedB = diag((omega*omega - N2(solver.zNative))/g)*solver.T;

            interiorRows = 2:(nEVP-1);
            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(A(1,:), solver.T(1,:), AbsTol=1e-12)
            testCase.verifyEqual(B(1,:), zeros(1,nEVP), AbsTol=1e-12)
            testCase.verifyEqual(A(end,:), solver.T(end,:), AbsTol=1e-12)
            testCase.verifyEqual(B(end,:), zeros(1,nEVP), AbsTol=1e-12)
        end

        function hydrostaticGModesUseZeroFrequencyWaveEVP(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            hydrostaticEVP = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g);
            zeroFrequencyEVP = IMEigenvalueProblem.waveModesAtFrequency(omega=0, f0=f0, g=g);

            [hydrostaticA, hydrostaticB] = hydrostaticEVP.assemble(solver);
            [zeroFrequencyA, zeroFrequencyB] = zeroFrequencyEVP.assemble(solver);

            testCase.verifyEqual(hydrostaticEVP.parameters.problemType, "waveModesAtFrequency")
            testCase.verifyEqual(hydrostaticEVP.parameters.omega, 0, AbsTol=0)
            testCase.verifyEqual(hydrostaticA, zeroFrequencyA, AbsTol=0)
            testCase.verifyEqual(hydrostaticB, zeroFrequencyB, AbsTol=0)
        end

        function waveModeMetadataSeparatesInnerWeightsAndNormalizations(testCase)
            [~, ~, ~, f0, g] = testCase.profile();
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);

            testCase.verifyTrue(isfield(evp.components.G, "innerWeight"))
            testCase.verifyTrue(isa(evp.components.G.innerWeight, "function_handle"))
            testCase.verifyEqual(evp.components.G.surfaceWeight, 1, AbsTol=0)
            testCase.verifyEqual(evp.components.G.bottomWeight, 0, AbsTol=0)
            testCase.verifyFalse(isfield(evp.components.G, "gramWeight"))
            testCase.verifyFalse(isfield(evp.components.G, "normalizationWeights"))
            testCase.verifyFalse(isfield(evp.components.G, "normalizationBoundaryWeights"))
            testCase.verifyTrue(isfield(evp.normalizations, "unity"))
            testCase.verifyTrue(isfield(evp.normalizations, "kConstant"))
            testCase.verifyTrue(isa(evp.normalizations.unity, "function_handle"))
        end

        function wkbCoordinateAppliesJacobianToSecondDerivative(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverWKBSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            z = solver.zNative;
            N = sqrt(N2(z));
            Nz = gradient(sqrt(N2(solver.zReference)), solver.zReference);
            Nz = interp1(solver.zReference, Nz, z, "pchip");

            D2 = solver.physicalDerivativeMatrix(2);
            expectedD2 = diag(N.*N)*solver.Txx + diag(Nz)*solver.Tx;

            testCase.verifyLessThan(norm(D2 - expectedD2, "fro")/norm(expectedD2, "fro"), 1e-10)
        end

        function densityCoordinateAppliesJacobianToSecondDerivative(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverDensitySpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            z = solver.zNative;
            N2z = gradient(N2(solver.zReference), solver.zReference);
            N2z = interp1(solver.zReference, N2z, z, "pchip");

            D2 = solver.physicalDerivativeMatrix(2);
            expectedD2 = diag(N2(z).*N2(z))*solver.Txx + diag(N2z)*solver.Tx;

            testCase.verifyLessThan(norm(D2 - expectedD2, "fro")/norm(expectedD2, "fro"), 1e-10)
        end

        function spectralUnityNormalizationUsesChebyshevInnerProduct(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            expectedFactors = testCase.expectedSpectralFactors(solver, basisSet, evp, ones(size(solver.zNative)));

            testCase.verifyEqual(basisSet.normalizationFactors(Normalization.unity), expectedFactors, RelTol=1e-11)
        end

        function wkbUnityNormalizationUsesNativeJacobian(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverWKBSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            expectedFactors = testCase.expectedSpectralFactors(solver, basisSet, evp, sqrt(N2(solver.zNative)));

            testCase.verifyEqual(basisSet.normalizationFactors(Normalization.unity), expectedFactors, RelTol=1e-11)
        end

        function densityUnityNormalizationUsesNativeJacobian(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverDensitySpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            expectedFactors = testCase.expectedSpectralFactors(solver, basisSet, evp, N2(solver.zNative));

            testCase.verifyEqual(basisSet.normalizationFactors(Normalization.unity), expectedFactors, RelTol=1e-11)
        end

        function solveEVPReturnsBasisWithChangeableNormalization(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);
            z = linspace(zDomain(1), zDomain(2), 24).';

            basisSet.normalization = Normalization.kConstant;
            Gk = basisSet.evaluate("G", z);
            basisSet.normalization = Normalization.wMax;
            Gw = basisSet.evaluate("G", z);

            testCase.verifySize(Gk, [length(z) 4])
            testCase.verifySize(Gw, [length(z) 4])
            testCase.verifyGreaterThan(norm(Gk - Gw, "fro"), 0)
        end

        function finiteDifferenceSolverUsesSuppliedGrid(testCase)
            [N2, zDomain, ~, f0, g] = testCase.profile();
            z = linspace(zDomain(1), zDomain(2), 18).';
            solver = IMSolverFiniteDifference(z=z, N2=N2, f0=f0, g=g);

            testCase.verifyEqual(solver.zNative, sort(z, "descend"))
            testCase.verifyEqual(solver.physicalDerivativeMatrix(0), eye(length(z)), AbsTol=0)
            testCase.verifyEqual(solver.boundaryIndex("surface"), 1)
            testCase.verifyEqual(solver.boundaryIndex("bottom"), length(z))
            testCase.verifyTrue(isa(solver, "IMSolver"))
            basisSet = solver.solveEVP(IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g), nModes=3);
            testCase.verifyClass(basisSet, "IMBasisSet")
            testCase.verifyNumElements(basisSet.h, 3)
        end

        function partialDepthPEIndexPolicyCountsBoundaryModes(testCase)
            positivePolicy = IMEigenvalueProblem.partialDepthPEIndexPolicy(boundarySign="positive");
            positiveIndex = positivePolicy.classify([1; 2; 3], struct());

            negativePolicy = IMEigenvalueProblem.partialDepthPEIndexPolicy(boundarySign="negative");
            negativeIndex = negativePolicy.classify([-2; -1; 3], struct());

            testCase.verifyEqual(positiveIndex.expectedNegativeCount, 0)
            testCase.verifyEqual(positiveIndex.negativeCount, 0)
            testCase.verifyEqual(negativeIndex.expectedNegativeCount, 2)
            testCase.verifyEqual(negativeIndex.negativeCount, 2)
        end

        function indexPolicyErrorsWhenObservedIndexDisagrees(testCase)
            policy = IMEigenvalueProblem.partialDepthPEIndexPolicy(boundarySign="negative");

            testCase.verifyError(@() policy.classify([-1; 2; 3], struct()), ...
                "IMIndexPolicy:IndexMismatch")
        end

        function hydrostaticFModeEVPDeclaresBarotropicZeroIndex(testCase)
            [~, ~, ~, ~, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);

            index = evp.indexPolicy.classify([0; 1; 2; 3], struct());

            testCase.verifyEqual(index.expectedZeroCount, 1)
            testCase.verifyEqual(index.zeroCount, 1)
            testCase.verifyEqual(index.positiveCount, 3)
        end

        function fullDepthPartialGramMatchesGramMatrix(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            fullGram = basisSet.gramMatrix("G");
            partialGram = basisSet.partialGramMatrix("G", zDomain(1), zDomain(2));

            testCase.verifyEqual(partialGram, fullGram, AbsTol=0)
            testCase.verifyEqual(partialGram, partialGram.', AbsTol=1e-12)
        end

        function unityNormalizationMakesSpectralGramUnitDiagonal(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            basisSet.normalization = Normalization.unity;
            gram = basisSet.gramMatrix("G");

            testCase.verifyEqual(diag(gram).', ones(1,4), AbsTol=1e-10)
        end

        function freeSurfaceInnerProductAddsSurfaceWeightOutsideIntegral(testCase)
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            zDomain = [-1300 0];
            nEVP = 48;
            g = 9.81;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=0, g=g);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0, g=g, upperBoundary="freeSurface");
            basisSet = solver.solveEVP(evp, nModes=3);
            z = solver.innerProductGrid(zDomain);
            rawG = solver.T*basisSet.nativeModes;
            weight = IMOperator.evaluateCoefficient(evp.components.G.innerWeight, z, testCase.evpContext(solver, evp));
            interior = solver.integrateInnerProduct(z, weight.*rawG(:,1).*rawG(:,1), zDomain);
            expectedFactor = sqrt(abs(interior + evp.components.G.surfaceWeight*rawG(1,1)*rawG(1,1)));

            factors = basisSet.normalizationFactors(Normalization.unity);
            basisSet.normalization = Normalization.unity;
            fullGram = basisSet.gramMatrix("G");
            zCut = zDomain(2) - 1;
            partialGram = basisSet.partialGramMatrix("G", zDomain(1), zCut);
            GPartial = basisSet.evaluate("G", z);
            expectedPartial = solver.integrateInnerProduct(z, weight.*GPartial(:,1).*GPartial(:,1), [zDomain(1) zCut]);

            testCase.verifyEqual(factors(1), expectedFactor, RelTol=1e-10)
            testCase.verifyEqual(fullGram(1,1), 1, AbsTol=1e-10)
            testCase.verifyEqual(partialGram(1,1), expectedPartial, AbsTol=1e-10)
            testCase.verifyGreaterThan(fullGram(1,1), partialGram(1,1))
        end

        function finiteDifferenceUnityNormalizationUsesBoundedTrapz(testCase)
            [N2, zDomain, ~, f0, g] = testCase.profile();
            solver = IMSolverFiniteDifference(z=linspace(zDomain(1), zDomain(2), 35).', N2=N2, f0=f0, g=g);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g, upperBoundary="freeSurface");
            basisSet = solver.solveEVP(evp, nModes=3);
            z = solver.innerProductGrid(zDomain);
            rawG = solver.evaluateNativeModes(basisSet.nativeModes, z);
            weight = IMOperator.evaluateCoefficient(evp.components.G.innerWeight, z, testCase.evpContext(solver, evp));
            normValue = trapz(z, weight.*rawG(:,1).*rawG(:,1));
            surfaceValue = solver.evaluateNativeModes(basisSet.nativeModes(:,1), zDomain(2));
            expectedFactor = sqrt(abs(normValue + evp.components.G.surfaceWeight*surfaceValue*surfaceValue));

            factors = basisSet.normalizationFactors(Normalization.unity);

            testCase.verifyEqual(factors(1), expectedFactor, RelTol=1e-10)
        end

        function constantStratificationFactoryReturnsAnalyticalBasisSet(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            nModes = 6;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0);

            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);

            testCase.verifyClass(basisSet, "IMBasisSetConstantStratification")
            testCase.verifyNumElements(basisSet.h, nModes)
            testCase.verifyEqual(basisSet.eigenvalues, 1./basisSet.h, RelTol=1e-12)
        end

        function constantStratificationInfersParametersFromEVP(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            g = 10;
            k = 1e-4;
            zDomain = [-1300 0];
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=f0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=2);
            D = diff(zDomain);
            hExpected = (N0*N0 - f0*f0)/(g*(k*k + (pi/D)^2));

            testCase.verifyEqual(basisSet.f0, f0, AbsTol=0)
            testCase.verifyEqual(basisSet.g, g, AbsTol=0)
            testCase.verifyEqual(basisSet.h(1), hExpected, RelTol=1e-12)
        end

        function constantStratificationMatchesDirectFixedWavenumberSolution(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            z = linspace(zDomain(1), zDomain(2), 65).';
            nModes = 6;
            k = 1e-4;
            g = 9.81;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes, g=g);

            direct = InternalModesConstantStratification(N0=N0, zIn=zDomain, zOut=z, latitude=0, nModes=nModes, g=g);
            direct.normalization = Normalization.kConstant;
            [FExpected, GExpected, hExpected] = direct.modesAtWavenumber(k);

            testCase.verifyEqual(basisSet.h, hExpected, RelTol=1e-12)
            testCase.verifyEqual(basisSet.evaluate("G", z), GExpected, AbsTol=1e-10)
            testCase.verifyEqual(basisSet.evaluate("F", z), FExpected, AbsTol=1e-10)
        end

        function constantStratificationMatchesDirectFixedFrequencySolution(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            z = linspace(zDomain(1), zDomain(2), 65).';
            nModes = 6;
            omega = 0.8*N0;
            g = 9.81;
            evp = IMEigenvalueProblem.waveModesAtFrequency(omega=omega, f0=0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes, g=g);

            direct = InternalModesConstantStratification(N0=N0, zIn=zDomain, zOut=z, latitude=0, nModes=nModes, g=g);
            direct.normalization = Normalization.kConstant;
            [FExpected, GExpected, hExpected] = direct.modesAtFrequency(omega);

            testCase.verifyEqual(basisSet.h, hExpected, RelTol=1e-12)
            testCase.verifyEqual(basisSet.evaluate("G", z), GExpected, AbsTol=1e-10)
            testCase.verifyEqual(basisSet.evaluate("F", z), FExpected, AbsTol=1e-10)
        end

        function constantStratificationFreeSurfaceMatchesDirectSolution(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            z = linspace(zDomain(1), zDomain(2), 65).';
            nModes = 6;
            k = 1e-4;
            g = 9.81;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=0, g=g, upperBoundary="freeSurface");
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes, g=g);

            direct = InternalModesConstantStratification(N0=N0, zIn=zDomain, zOut=z, latitude=0, nModes=nModes, g=g);
            direct.upperBoundary = UpperBoundary.freeSurface;
            direct.normalization = Normalization.kConstant;
            [FExpected, GExpected, hExpected] = direct.modesAtWavenumber(k);

            testCase.verifyEqual(basisSet.h, hExpected, RelTol=1e-10)
            testCase.verifyEqual(basisSet.evaluate("G", z), GExpected, AbsTol=1e-8)
            testCase.verifyEqual(basisSet.evaluate("F", z), FExpected, AbsTol=1e-8)
        end

        function constantStratificationNormalizationCanChange(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            z = linspace(zDomain(1), zDomain(2), 65).';
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=4);

            basisSet.normalization = Normalization.kConstant;
            Gk = basisSet.evaluate("G", z);
            basisSet.normalization = Normalization.wMax;
            Gw = basisSet.evaluate("G", z);

            testCase.verifyGreaterThan(norm(Gk - Gw, "fro"), 0)
        end

        function constantStratificationGramUsesEVPInnerWeight(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            nModes = 4;
            g = 9.81;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes, g=g);
            z = linspace(zDomain(1), zDomain(2), 256).';
            G = basisSet.evaluate("G", z);
            weight = (N0*N0/g)*ones(size(z));
            expectedGram = zeros(nModes,nModes);
            for iMode = 1:nModes
                for jMode = iMode:nModes
                    value = trapz(z, weight.*G(:,iMode).*G(:,jMode));
                    expectedGram(iMode,jMode) = value;
                    expectedGram(jMode,iMode) = value;
                end
            end

            testCase.verifyEqual(basisSet.gramMatrix("G"), expectedGram, AbsTol=1e-12)
        end

        function unityNormalizationMakesConstantStratificationGramUnitDiagonal(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            nModes = 4;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);

            basisSet.normalization = Normalization.unity;
            gram = basisSet.gramMatrix("G");

            testCase.verifyEqual(diag(gram).', ones(1,nModes), AbsTol=1e-12)
            testCase.verifyLessThan(norm(gram - diag(diag(gram)), "fro"), 1e-12)
        end

        function deferredAnalyticalFactoriesThrowExplicitErrors(testCase)
            testCase.verifyError(@() IMBasisSet.exponentialStratification(), ...
                "IMBasisSet:AnalyticalBasisNotImplemented")
            testCase.verifyError(@() IMBasisSet.wkbApproximation(), ...
                "IMBasisSet:AnalyticalBasisNotImplemented")
        end
    end

    methods (Access = private)
        function [N2, zDomain, nEVP, f0, g] = profile(~)
            N0 = 5.2e-3;
            b = 1300;
            N2 = @(z) N0*N0*exp(2*z/b);
            zDomain = [-4000 0];
            nEVP = 32;
            f0 = 1e-4;
            g = 9.81;
        end

        function factors = expectedSpectralFactors(testCase, solver, basisSet, evp, q)
            z = solver.zNative;
            rawG = solver.T*basisSet.nativeModes;
            weight = IMOperator.evaluateCoefficient(evp.components.G.innerWeight, z, testCase.evpContext(solver, evp));
            weights = testCase.chebyshevIntegrationWeights(solver);
            factors = zeros(1,size(rawG,2));
            for iMode = 1:size(rawG,2)
                integrandCheb = InternalModesSpectral.fct(weight.*rawG(:,iMode).*rawG(:,iMode)./q(:));
                normValue = sum(weights.*integrandCheb);
                normValue = normValue + evp.components.G.surfaceWeight*rawG(1,iMode)*rawG(1,iMode);
                normValue = normValue + evp.components.G.bottomWeight*rawG(end,iMode)*rawG(end,iMode);
                factors(iMode) = sqrt(abs(normValue));
            end
        end

        function context = evpContext(~, solver, evp)
            context = solver.context();
            names = fieldnames(evp.parameters);
            for iName = 1:length(names)
                context.(names{iName}) = evp.parameters.(names{iName});
            end
        end

        function weights = chebyshevIntegrationWeights(~, solver)
            np = (0:(solver.nEVP-1)).';
            weights = -(1+(-1).^np)./(np.*np - 1);
            weights(2) = 0;
            weights = (max(solver.xNative) - min(solver.xNative))*weights/2;
        end
    end
end
