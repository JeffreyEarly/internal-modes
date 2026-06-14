classdef IMEigenvalueProblemEndpointCaseTests < matlab.unittest.TestCase

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
        function endpointAlgebraUsesYassinIndexing(testCase)
            surface = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            bottom = IMBoundaryCondition(a=0, b=1, c=1, d=0);

            testCase.verifyTrue(surface.isEigenvalueDependent())
            testCase.verifyEqual(surface.determinant("surface"), 1, AbsTol=0)
            testCase.verifyEqual(surface.metricWeight("surface"), 1, AbsTol=0)
            testCase.verifyEqual(bottom.determinant("bottom"), -1, AbsTol=0)
            testCase.verifyEqual(bottom.metricWeight("bottom"), -1, AbsTol=0)
        end

        function tinyEigenvalueSideCoefficientsAreInactive(testCase)
            boundary = IMBoundaryCondition(a=1, b=0, c=eps, d=0);
            evp = IMEigenvalueProblem(surfaceBoundary=boundary);

            testCase.verifyFalse(boundary.isEigenvalueDependent())
            testCase.verifyEmpty(evp.endpointWeights("surface"))
        end

        function inactiveRobinEnergyUsesYassinEndpointSigns(testCase)
            boundary = IMBoundaryCondition(a=2, b=3);

            testCase.verifyEqual(boundary.robinEnergyCoefficient("surface"), -2/3, AbsTol=0)
            testCase.verifyEqual(boundary.robinEnergyCoefficient("bottom"), 2/3, AbsTol=0)
        end

        function case1PositiveMetricAndNonnegativeQSuppressesNegativeCandidates(testCase)
            solver = testCase.canonicalSolver();
            evp = IMEigenvalueProblem(p=1, q=0, r=1);
            [A, ~] = evp.assemble(solver);

            info = evp.definitenessInfo(solver);
            bounds = evp.negativeEigenvalueBounds(solver, A);
            selection = evp.selectModes([-100; 1; 2; 3], 2, solver, A);

            testCase.verifyTrue(info.metricPositive)
            testCase.verifyTrue(info.qNonnegativeCertified)
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 0)
            testCase.verifyEqual(selection.sortIndex(:), [2; 3])
            testCase.verifyEqual(selection.modeNumber, [1 2])
        end

        function case2MetricIndexWithZeroAbsentHasExactNegativeCount(testCase)
            solver = testCase.canonicalSolver();
            bottom = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMEigenvalueProblem(p=1, q=1, r=1, bottomBoundary=bottom);
            [A, ~] = evp.assemble(solver);

            info = evp.definitenessInfo(solver);
            bounds = evp.negativeEigenvalueBounds(solver, A);

            testCase.verifyEqual(info.metricIndex, 1)
            testCase.verifyTrue(info.qNonnegativeCertified)
            testCase.verifyEqual(bounds.zeroEigenvalueStatus, "absent")
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 1)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 1)
        end

        function case2WithoutLeftMatrixKeepsMetricIndexAsUpperBound(testCase)
            solver = testCase.canonicalSolver();
            bottom = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMEigenvalueProblem(p=1, q=1, r=1, bottomBoundary=bottom);

            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyEqual(bounds.zeroEigenvalueStatus, "unchecked")
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 1)
        end

        function case3InactiveRobinDirectionsSetFiniteSearchCap(testCase)
            solver = testCase.canonicalSolver();
            surface = IMBoundaryCondition(a=1, b=1);
            bottom = IMBoundaryCondition(a=-1, b=1);
            evp = IMEigenvalueProblem(p=1, q=0, r=1, surfaceBoundary=surface, bottomBoundary=bottom);

            info = evp.definitenessInfo(solver);
            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyTrue(info.metricPositive)
            testCase.verifyEqual(info.endpointNumeratorNegativeDirections, 2)
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 2)
        end

        function case3ActiveEndpointInertiaSetsFiniteSearchCap(testCase)
            solver = testCase.canonicalSolver();
            surface = IMBoundaryCondition(a=0, b=1, c=1, d=1);
            evp = IMEigenvalueProblem(p=1, q=0, r=1, surfaceBoundary=surface);

            info = evp.definitenessInfo(solver);
            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyTrue(info.metricPositive)
            testCase.verifyEqual(info.endpointNumeratorNegativeDirections, 1)
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 1)
        end

        function case3GeostrophicWKBRetainsCertifiedNegativeMode(testCase)
            [solver, evp] = testCase.geostrophicWKBProblem();
            [A, B] = evp.assemble(solver);
            eigenvalues = eig(A, B);
            valid = isfinite(real(eigenvalues)) & isfinite(imag(eigenvalues)) ...
                & abs(imag(eigenvalues)) < 1e-8*max(1,abs(real(eigenvalues)));
            finiteRealEigenvalues = real(eigenvalues(valid));

            info = evp.definitenessInfo(solver);
            bounds = evp.negativeEigenvalueBounds(solver, A);
            basisSet = solver.solveEVP(evp, nModes=4);

            testCase.verifyGreaterThanOrEqual(nnz(finiteRealEigenvalues < 0), 1)
            testCase.verifyEqual(info.endpointNumeratorNegativeDirections, 2)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 2)
            testCase.verifyLessThan(basisSet.eigenvalues(1), 0)
            testCase.verifyEqual(basisSet.modeNumber(1), -1)
        end

        function case4FailedCoefficientSignsReturnUnknown(testCase)
            solver = testCase.canonicalSolver();
            evp = IMEigenvalueProblem(p=-1, q=0, r=1);

            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, "unknown")
            testCase.verifyTrue(contains(bounds.reason, "fail"))
        end

        function case4NonfiniteSamplesReturnUnknown(testCase)
            solver = testCase.canonicalSolver();
            evp = IMEigenvalueProblem(p=@(z,~) NaN(size(z)), q=0, r=1);

            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, "unknown")
            testCase.verifyTrue(contains(bounds.reason, "nonfinite"))
        end

        function case4DegenerateActiveEndpointReturnsUnknown(testCase)
            solver = testCase.canonicalSolver();
            surface = IMBoundaryCondition(a=0, b=1, c=0, d=1);
            evp = IMEigenvalueProblem(p=1, q=0, r=1, surfaceBoundary=surface);

            info = evp.definitenessInfo(solver);
            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyTrue(info.hasDegenerateEndpointMetric)
            testCase.verifyEmpty(evp.endpointWeights("surface"))
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, "unknown")
            testCase.verifyTrue(contains(bounds.reason, "degenerate"))
        end
    end

    methods (Static, Access = private)
        function solver = canonicalSolver()
            N2 = @(z) ones(size(z));
            solver = IMSolverSpectral(N2=N2, zDomain=[-1 0], nEVP=24);
        end

        function [solver, evp] = geostrophicWKBProblem()
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            g = 9.81;
            g0 = -N0*N0*b;
            gd = -N0*N0*b/10;
            zDomain = [-D 0];
            N2 = @(z) N0*N0*exp(2*z/b);

            p = @(z,ctx) 1 ./ ctx.N2(z);
            q = @(z,~) zeros(size(z));
            r = @(z,ctx) ones(size(z))/ctx.g;
            surfaceBoundary = IMBoundaryCondition(a=-(1/g + 1/g0), b=1);
            bottomBoundary = IMBoundaryCondition(a=1/gd, b=1);
            normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);

            evp = IMInternalModes(name="unforced-APV-modes", formulation="F", ...
                p=p, q=q, r=r, g=g, normalizations=normalizations, ...
                defaultNormalization=Normalization.unity, surfaceBoundary=surfaceBoundary, ...
                bottomBoundary=bottomBoundary);
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=128, coordinateKind="wkb");
        end
    end
end
