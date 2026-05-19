classdef InternalModesEVPRefactorTests < matlab.unittest.TestCase

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
        function spectralHydrostaticGAssemblyMatchesPhysicalStrongForm(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            k = 1e-4;
            solver = InternalModesSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = InternalModesEVP.hydrostaticGModes(k=k, f0=f0, g=g);

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

        function wkbCoordinateAppliesJacobianToSecondDerivative(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = InternalModesSolverWKBSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
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
            solver = InternalModesSolverDensitySpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            z = solver.zNative;
            N2z = gradient(N2(solver.zReference), solver.zReference);
            N2z = interp1(solver.zReference, N2z, z, "pchip");

            D2 = solver.physicalDerivativeMatrix(2);
            expectedD2 = diag(N2(z).*N2(z))*solver.Txx + diag(N2z)*solver.Tx;

            testCase.verifyLessThan(norm(D2 - expectedD2, "fro")/norm(expectedD2, "fro"), 1e-10)
        end

        function solveEVPReturnsBasisWithChangeableNormalization(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = InternalModesSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = InternalModesEVP.hydrostaticGModes(k=1e-4, f0=f0, g=g);
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
            solver = InternalModesSolverFiniteDifference(z=z, N2=N2, f0=f0, g=g);

            testCase.verifyEqual(solver.zNative, sort(z, "descend"))
            testCase.verifyEqual(solver.physicalDerivativeMatrix(0), eye(length(z)), AbsTol=0)
            testCase.verifyEqual(solver.boundaryIndex("surface"), 1)
            testCase.verifyEqual(solver.boundaryIndex("bottom"), length(z))
        end

        function partialDepthPEIndexPolicyCountsBoundaryModes(testCase)
            positivePolicy = InternalModesEVP.partialDepthPEIndexPolicy(boundarySign="positive");
            positiveIndex = positivePolicy.classify([1; 2; 3], struct());

            negativePolicy = InternalModesEVP.partialDepthPEIndexPolicy(boundarySign="negative");
            negativeIndex = negativePolicy.classify([-2; -1; 3], struct());

            testCase.verifyEqual(positiveIndex.expectedNegativeCount, 0)
            testCase.verifyEqual(positiveIndex.negativeCount, 0)
            testCase.verifyEqual(negativeIndex.expectedNegativeCount, 2)
            testCase.verifyEqual(negativeIndex.negativeCount, 2)
        end

        function indexPolicyErrorsWhenObservedIndexDisagrees(testCase)
            policy = InternalModesEVP.partialDepthPEIndexPolicy(boundarySign="negative");

            testCase.verifyError(@() policy.classify([-1; 2; 3], struct()), ...
                "InternalModesIndexPolicy:IndexMismatch")
        end

        function hydrostaticFModeEVPDeclaresBarotropicZeroIndex(testCase)
            [~, ~, ~, ~, g] = testCase.profile();
            evp = InternalModesEVP.hydrostaticFModes(g=g);

            index = evp.indexPolicy.classify([0; 1; 2; 3], struct());

            testCase.verifyEqual(index.expectedZeroCount, 1)
            testCase.verifyEqual(index.zeroCount, 1)
            testCase.verifyEqual(index.positiveCount, 3)
        end

        function fullDepthPartialGramMatchesGramMatrix(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = InternalModesSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, f0=f0, g=g);
            evp = InternalModesEVP.hydrostaticGModes(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            fullGram = basisSet.gramMatrix("G");
            partialGram = basisSet.partialGramMatrix("G", zDomain(1), zDomain(2));

            testCase.verifyEqual(partialGram, fullGram, AbsTol=0)
            testCase.verifyEqual(partialGram, partialGram.', AbsTol=1e-12)
        end

        function analyticalPlaceholderThrowsUnsupportedOperation(testCase)
            basisSet = InternalModesBasisSet.constantStratification();

            testCase.verifyError(@() basisSet.evaluate("G", linspace(-1,0,8).'), ...
                "InternalModesBasisSet:UnsupportedOperation")
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
    end
end
