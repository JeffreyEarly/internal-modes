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
        function canonicalAssemblyMatchesStrongFormForConstantCoefficients(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem(name="constant", p=2, q=3, r=4);

            [A, B] = evp.assemble(solver);
            expectedA = -2*solver.physicalDerivativeMatrix(2) + 3*solver.physicalDerivativeMatrix(0);
            expectedB = 4*solver.physicalDerivativeMatrix(0);
            interiorRows = 2:(nEVP-1);

            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-11)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-11)
            testCase.verifyEqual(A(1,:), -solver.T(1,:), AbsTol=1e-11)
            testCase.verifyEqual(B(1,:), zeros(1,nEVP), AbsTol=1e-11)
            testCase.verifyEqual(A(end,:), -solver.T(end,:), AbsTol=1e-11)
            testCase.verifyEqual(B(end,:), zeros(1,nEVP), AbsTol=1e-11)
        end

        function canonicalAssemblyIncludesGridDerivativeOfVariableP(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            p = @(z,~) 1 + (z/1000).^2;
            evp = IMEigenvalueProblem(name="variableP", p=p, q=0, r=1);

            [A, B] = evp.assemble(solver);
            pValues = p(solver.zNative, struct());
            pzValues = solver.differentiateGridValues(pValues, 1);
            expectedA = -diag(pValues)*solver.physicalDerivativeMatrix(2) ...
                - diag(pzValues)*solver.physicalDerivativeMatrix(1);
            expectedB = solver.physicalDerivativeMatrix(0);
            interiorRows = 2:(nEVP-1);

            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-9)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-11)
        end

        function endpointRowsUseCanonicalBoundaryCoefficients(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            surface = IMBoundaryCondition(a=2, b=3);
            bottom = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMEigenvalueProblem(p=5, surfaceBoundary=surface, bottomBoundary=bottom);

            [A, B] = evp.assemble(solver);
            D0 = solver.physicalDerivativeMatrix(0);
            D1 = solver.physicalDerivativeMatrix(1);

            testCase.verifyEqual(A(1,:), -2*D0(1,:) + 15*D1(1,:), AbsTol=1e-11)
            testCase.verifyEqual(B(1,:), zeros(1,nEVP), AbsTol=1e-11)
            testCase.verifyEqual(A(end,:), 5*D1(end,:), AbsTol=1e-11)
            testCase.verifyEqual(B(end,:), D0(end,:), AbsTol=1e-11)
        end

        function boundaryConditionReportsMetricWeights(testCase)
            surface = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            bottom = IMBoundaryCondition(a=0, b=1, c=1, d=0);

            testCase.verifyTrue(surface.isEigenvalueDependent())
            testCase.verifyEqual(surface.determinant("surface"), -1, AbsTol=0)
            testCase.verifyEqual(surface.metricWeight("surface"), -1, AbsTol=0)
            testCase.verifyEqual(bottom.determinant("bottom"), 1, AbsTol=0)
            testCase.verifyEqual(bottom.metricWeight("bottom"), 1, AbsTol=0)
            testCase.verifyEqual(IMBoundaryCondition.dirichlet().a, 1, AbsTol=0)
            testCase.verifyEqual(IMBoundaryCondition.neumann().b, 1, AbsTol=0)
            testCase.verifyEqual(IMBoundaryCondition.robin(2,3).b, 3, AbsTol=0)
        end

        function definitenessDiagnosticsIgnoreGarbageNegativesWhenCertified(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem(p=1, q=0, r=1);
            [A, ~] = evp.assemble(solver);

            info = evp.definitenessInfo(solver);
            bounds = evp.negativeEigenvalueBounds(solver, A);
            selection = evp.selectModes([-100; 1; 2; 3], 2, solver, A);

            testCase.verifyTrue(info.metricPositive)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 0)
            testCase.verifyEqual(selection.sortIndex(:), [2; 3])
            testCase.verifyEqual(selection.modeNumber, [1 2])
        end

        function internalModeFactoriesAssembleCanonicalForms(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            k = 1e-4;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMInternalModes.waveModesAtWavenumber(k=k, f0=f0, g=g);

            [A, B] = evp.assemble(solver);
            expectedA = -solver.physicalDerivativeMatrix(2) + k*k*solver.physicalDerivativeMatrix(0);
            expectedB = diag((N2(solver.zNative) - f0*f0)/g)*solver.physicalDerivativeMatrix(0);
            interiorRows = 2:(nEVP-1);

            testCase.verifyClass(evp, "IMInternalModes")
            testCase.verifyEqual(evp.formulation, "G")
            testCase.verifyEqual(evp.metadata.k, k, AbsTol=0)
            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-11)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-11)
        end

        function solverReturnsInternalModesBasisWithFAndG(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes();

            basisSet = solver.solveEVP(evp, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 16).';

            testCase.verifyClass(basisSet, "IMInternalModesBasis")
            testCase.verifySize(basisSet.G(z), [16 3])
            testCase.verifySize(basisSet.F(z), [16 3])
            testCase.verifyEqual(size(basisSet.gramMatrix("G")), [3 3])
        end

        function spectralSolverSupportsCoordinateKinds(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            coordinateKinds = ["z", "wkb", "density"];
            for coordinateKind = coordinateKinds
                solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP, coordinateKind=coordinateKind);
                testCase.verifyEqual(solver.coordinateKind, coordinateKind)
                testCase.verifySize(solver.physicalDerivativeMatrix(1), [nEVP nEVP])
            end
        end

        function finiteDifferenceSolverSolvesCanonicalProblem(testCase)
            z = linspace(-1000, 0, 32).';
            solver = IMSolverFiniteDifference(z=z, N2=@(z) 1e-5*ones(size(z)));
            evp = IMEigenvalueProblem(p=1, q=0, r=1);

            basisSet = solver.solveEVP(evp, nModes=2);

            testCase.verifyClass(basisSet, "IMBasisSet")
            testCase.verifySize(basisSet.u(linspace(-1000,0,8).'), [8 2])
        end

        function noValidEigenvalueDiagnosticStillReportsMatrixStats(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem(name="degenerate", p=0, q=0, r=0);

            testCase.verifyError(@() solver.solveEVP(evp, nModes=2), "IMSolver:NoValidEigenvalues")
        end
    end

    methods (Static, Access = private)
        function [N2, zDomain, nEVP, f0, g] = profile()
            N0 = 5.2e-3;
            zDomain = [-1000 0];
            nEVP = 24;
            f0 = 1e-4;
            g = 9.81;
            N2 = @(z) N0*N0*(1 + 0.1*z/abs(zDomain(1)));
        end
    end
end
