classdef IMSurfaceGeostrophicModesTests < matlab.unittest.TestCase

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
        function spectralSurfaceModesMatchConstantSolution(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            zDomain = [-5000 0];
            k = [1e-4 2e-4];
            N2 = @(z) N0*N0*ones(size(z));
            problem = IMSurfaceGeostrophicModes.surfaceModesAtWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=k);
            solver = IMSolverSpectral(nEVP=96);
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, f0=f0);

            basisSet = solver.solveSurfaceGeostrophicModes(problem);
            exact = solution.sqgModesAtWavenumber(k, boundary="surface");
            z = linspace(zDomain(1), zDomain(2), 17).';

            testCase.verifyClass(basisSet, "IMSurfaceGeostrophicModesBasis")
            testCase.verifyEqual(basisSet.psi(z), exact.psi(z), AbsTol=1e-6)
        end

        function spectralBottomModesMatchConstantSolution(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            zDomain = [-5000 0];
            k = [1e-4 2e-4];
            N2 = @(z) N0*N0*ones(size(z));
            problem = IMSurfaceGeostrophicModes.bottomModesAtWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=k);
            solver = IMSolverSpectral(nEVP=96);
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, f0=f0);

            basisSet = solver.solveSurfaceGeostrophicModes(problem);
            exact = solution.sqgModesAtWavenumber(k, boundary="bottom");
            z = linspace(zDomain(1), zDomain(2), 17).';

            testCase.verifyEqual(basisSet.psi(z), exact.psi(z), AbsTol=1e-6)
        end

        function spectralModesMatchExponentialSolution(testCase)
            N0 = 5.2e-3;
            b = 1300;
            f0 = 1e-4;
            zDomain = [-5000 0];
            k = [1e-4 2e-4];
            N2 = @(z) N0*N0*exp(2*z/b);
            solver = IMSolverSpectral(nEVP=160);
            solution = IMExponentialStratificationSolution(N0=N0, b=b, zDomain=zDomain, f0=f0);
            z = linspace(zDomain(1), zDomain(2), 17).';

            surface = solver.solveSurfaceGeostrophicModes(IMSurfaceGeostrophicModes.surfaceModesAtWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=k));
            bottom = solver.solveSurfaceGeostrophicModes(IMSurfaceGeostrophicModes.bottomModesAtWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=k));

            testCase.verifyEqual(surface.psi(z), solution.sqgModesAtWavenumber(k, boundary="surface").psi(z), AbsTol=1e-4)
            testCase.verifyEqual(bottom.psi(z), solution.sqgModesAtWavenumber(k, boundary="bottom").psi(z), AbsTol=1e-4)
        end

        function vectorWavenumbersSetOutputColumns(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            zDomain = [-1000 0];
            k = [1e-4 2e-4 3e-4];
            N2 = @(z) N0*N0*ones(size(z));
            problem = IMSurfaceGeostrophicModes(N2=N2, zDomain=zDomain, f0=f0, k=k, boundary="surface");
            solver = IMSolverSpectral(nEVP=64);
            basisSet = solver.solveSurfaceGeostrophicModes(problem);
            z = linspace(zDomain(1), zDomain(2), 11).';

            testCase.verifyEqual(basisSet.k, k)
            testCase.verifySize(basisSet.psi(z), [length(z) length(k)])
            testCase.verifySize(basisSet.psiz(z), [length(z) length(k)])
        end

        function boundaryForcingIsAppliedOnSolverGrid(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            zDomain = [-5000 0];
            k = [1e-4 2e-4];
            N2 = @(z) N0*N0*ones(size(z));
            solver = IMSolverSpectral(nEVP=96);

            surface = solver.solveSurfaceGeostrophicModes(IMSurfaceGeostrophicModes.surfaceModesAtWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=k));
            bottom = solver.solveSurfaceGeostrophicModes(IMSurfaceGeostrophicModes.bottomModesAtWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=k));

            testCase.verifyEqual(f0*surface.psiz(zDomain(2)), ones(1,length(k)), AbsTol=1e-9)
            testCase.verifyEqual(f0*surface.psiz(zDomain(1)), zeros(1,length(k)), AbsTol=1e-9)
            testCase.verifyEqual(f0*bottom.psiz(zDomain(2)), zeros(1,length(k)), AbsTol=1e-9)
            testCase.verifyEqual(f0*bottom.psiz(zDomain(1)), ones(1,length(k)), AbsTol=1e-9)
        end

        function finiteDifferenceSolvesAndValidatesDomain(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            zDomain = [-1000 0];
            N2 = @(z) N0*N0*ones(size(z));
            problem = IMSurfaceGeostrophicModes.surfaceModesAtWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=1e-4);
            solver = IMSolverFiniteDifference(z=linspace(zDomain(1), zDomain(2), 65).');

            basisSet = solver.solveSurfaceGeostrophicModes(problem);

            testCase.verifyClass(basisSet, "IMSurfaceGeostrophicModesBasis")
            testCase.verifySize(basisSet.psi(linspace(zDomain(1), zDomain(2), 8).'), [8 1])
            mismatchedSolver = IMSolverFiniteDifference(z=linspace(-900, 0, 65).');
            testCase.verifyError(@() mismatchedSolver.solveSurfaceGeostrophicModes(problem), "IMSolverFiniteDifference:DomainMismatch")
        end

        function zeroCoriolisIsRejected(testCase)
            N2 = @(z) ones(size(z));

            testCase.verifyError(@() IMSurfaceGeostrophicModes.surfaceModesAtWavenumber(N2=N2, zDomain=[-1 0], f0=0, k=1), "IMSurfaceGeostrophicModes:InvalidCoriolis")
        end
    end
end
