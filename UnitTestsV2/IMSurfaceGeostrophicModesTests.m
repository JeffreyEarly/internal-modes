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
        function constructorSupportsEndpointCombinations(testCase)
            N2 = @(z) ones(size(z));

            surface = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=[-1 0], f0=1, k=2, g0=3);
            bottom = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=[-1 0], f0=1, k=2, gd=4);
            both = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=[-1 0], f0=1, k=2, g0=3, gd=4, surfaceAnomaly="noFreeSurface");

            testCase.verifyEqual(surface.modesPerWavenumber(), 1)
            testCase.verifyEqual(bottom.modesPerWavenumber(), 1)
            testCase.verifyEqual(both.modesPerWavenumber(), 2)
            testCase.verifyEqual(surface.g0, 3)
            testCase.verifyTrue(isinf(surface.gd))
            testCase.verifyEqual(bottom.gd, 4)
            testCase.verifyEqual(both.surfaceAnomaly, "noFreeSurface")
        end

        function constructorRejectsInvalidEndpointChoices(testCase)
            N2 = @(z) ones(size(z));

            testCase.verifyError(@() IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=[-1 0], f0=1, k=2), "IMSurfaceGeostrophicModes:NoBoundaryAnomaly")
            testCase.verifyError(@() IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=[-1 0], f0=1, k=2, g0=0), "IMSurfaceGeostrophicModes:InvalidEndpointWeight")
            testCase.verifyError(@() IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=[-1 0], f0=1, k=2, gd=0), "IMSurfaceGeostrophicModes:InvalidEndpointWeight")
            testCase.verifyError(@() IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=[-1 0], f0=0, k=2, g0=3), "IMSurfaceGeostrophicModes:InvalidCoriolis")
            testCase.verifyTrue(testCase.throwsAny(@() IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=[-1 0], f0=1, k=2, g0=3, surfaceAnomaly="bad")))
        end

        function surfaceOnlySolveReportsProjectedFields(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            g = 9.81;
            g0 = 0.04;
            zDomain = [-1000 0];
            k = [1e-4 2e-4];
            N2 = @(z) N0*N0*ones(size(z));
            problem = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=zDomain, f0=f0, g=g, k=k, g0=g0);
            solver = IMSolverSpectral(nEVP=96);

            basisSet = solver.solveSurfaceGeostrophicModes(problem);
            z = linspace(zDomain(1), zDomain(2), 11).';

            testCase.verifyClass(basisSet, "IMSurfaceGeostrophicModesBasis")
            testCase.verifyEqual(basisSet.k, k)
            testCase.verifyEqual(basisSet.modeNumber, [1 1])
            testCase.verifySize(basisSet.F(z), [length(z) length(k)])
            testCase.verifySize(basisSet.G(z), [length(z) length(k)])
            testCase.verifyEqual(basisSet.h, 2*basisSet.k.^2.*basisSet.energyEigenvalues, RelTol=1e-12)

            etaSurface = f0/g*(basisSet.G(zDomain(2)) - basisSet.F(zDomain(2)));
            etaBottom = f0/g*basisSet.G(zDomain(1));
            expectedSurface = basisSet.mixingCoefficients(1,:)*g0/(N0*N0);
            testCase.verifyEqual(etaSurface, expectedSurface, RelTol=1e-8, AbsTol=1e-9)
            testCase.verifyEqual(etaBottom, zeros(size(etaBottom)), AbsTol=1e-9)
        end

        function noFreeSurfaceSurfaceAnomalyOmitsFValue(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            g = 9.81;
            g0 = 0.04;
            zDomain = [-1000 0];
            N2 = @(z) N0*N0*ones(size(z));
            problem = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=zDomain, f0=f0, g=g, k=1e-4, g0=g0, surfaceAnomaly="noFreeSurface");
            basisSet = IMSolverSpectral(nEVP=96).solveSurfaceGeostrophicModes(problem);

            etaSurface = f0/g*basisSet.G(zDomain(2));
            expectedSurface = basisSet.mixingCoefficients(1,:)*g0/(N0*N0);
            testCase.verifyEqual(etaSurface, expectedSurface, RelTol=1e-8, AbsTol=1e-9)
        end

        function exponentialSurfaceSolveAvoidsNearlySingularWarning(testCase)
            N0 = 5.2e-3;
            b = 1300;
            f0 = 1e-4;
            g0 = -N0*N0*b;
            zDomain = [-4000 0];
            L = 1e6;
            k = 2*pi/L;
            N2 = @(z) N0*N0*exp(2*z/b);
            problem = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=k, g0=g0);
            solver = IMSolverSpectral(nEVP=256);

            warningState = warning("error", "MATLAB:nearlySingularMatrix");
            cleanup = onCleanup(@() warning(warningState));
            basisSet = solver.solveSurfaceGeostrophicModes(problem);
            z = linspace(zDomain(1), zDomain(2), 64).';
            F = basisSet.F(z);

            testCase.verifyTrue(all(isfinite(F(:))))
            testCase.verifyTrue(all(isfinite(basisSet.h)))
        end

        function twoBoundarySolveDiagonalizesBoundaryEnergy(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            g = 9.81;
            g0 = 0.04;
            gd = 0.02;
            zDomain = [-1000 0];
            N2 = @(z) N0*N0*ones(size(z));
            problem = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=zDomain, f0=f0, g=g, k=1e-4, g0=g0, gd=gd);

            basisSet = IMSolverSpectral(nEVP=96).solveSurfaceGeostrophicModes(problem);
            etaSurface = f0/g*(basisSet.G(zDomain(2)) - basisSet.F(zDomain(2)));
            etaBottom = f0/g*basisSet.G(zDomain(1));
            expectedSurface = basisSet.mixingCoefficients(1,:)*g0/(N0*N0);
            expectedBottom = basisSet.mixingCoefficients(2,:)*gd/(N0*N0);
            energyMatrix = -f0*(basisSet.F(zDomain(2)).'*etaSurface) + f0*(basisSet.F(zDomain(1)).'*etaBottom) + g0*(etaSurface.'*etaSurface) + gd*(etaBottom.'*etaBottom);
            energyMatrix = 0.5*(energyMatrix + energyMatrix.');

            testCase.verifyEqual(etaSurface, expectedSurface, RelTol=1e-8, AbsTol=1e-9)
            testCase.verifyEqual(etaBottom, expectedBottom, RelTol=1e-8, AbsTol=1e-9)
            testCase.verifyEqual(energyMatrix, diag(basisSet.energyEigenvalues), AbsTol=1e-9)
            testCase.verifyEqual(basisSet.modeNumber, [1 2])
            testCase.verifyEqual(basisSet.h, 2*basisSet.k.^2.*basisSet.energyEigenvalues, RelTol=1e-12)
        end

        function finiteDifferenceSolvesAndValidatesDomain(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            zDomain = [-1000 0];
            N2 = @(z) N0*N0*ones(size(z));
            problem = IMSurfaceGeostrophicModes.atWavenumber(N2=N2, zDomain=zDomain, f0=f0, k=1e-4, gd=0.02);
            solver = IMSolverFiniteDifference(z=linspace(zDomain(1), zDomain(2), 65).');

            basisSet = solver.solveSurfaceGeostrophicModes(problem);

            testCase.verifyClass(basisSet, "IMSurfaceGeostrophicModesBasis")
            testCase.verifySize(basisSet.F(linspace(zDomain(1), zDomain(2), 8).'), [8 1])
            mismatchedSolver = IMSolverFiniteDifference(z=linspace(-900, 0, 65).');
            testCase.verifyError(@() mismatchedSolver.solveSurfaceGeostrophicModes(problem), "IMSolverFiniteDifference:DomainMismatch")
        end
    end

    methods (Access = private)
        function didThrow = throwsAny(~, functionHandle)
            didThrow = false;
            try
                functionHandle();
            catch
                didThrow = true;
            end
        end
    end
end
