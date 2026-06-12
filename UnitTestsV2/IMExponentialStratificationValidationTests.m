classdef IMExponentialStratificationValidationTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryAndOraclePaths(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            parentRoot = fileparts(repoRoot);
            oracleRoot = fullfile(parentRoot, "cheb" + "fun");
            testCase.originalPath = path;
            addpath(repoRoot);
            if exist(oracleRoot, "dir")
                addpath(oracleRoot);
            end
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function factoryReturnsExponentialBasisSet(testCase)
            [N0, b, zDomain, ~, nModes, f0, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g);
            basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes);

            testCase.verifyClass(basisSet, "IMBasisSetExponentialStratification")
            testCase.verifyEqual(basisSet.N0, N0)
            testCase.verifyEqual(basisSet.b, b)
            testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
            testCase.verifyEqual(basisSet.N2(zDomain(2)), N0*N0, RelTol=1e-14)
            testCase.verifyGreaterThan(basisSet.F(zDomain(2)), zeros(1,nModes))
        end

        function fixedWavenumberGModesMatchV1Oracle(testCase)
            [N0, b, zDomain, z, nModes, f0, g] = testCase.profile();
            k = 1e-4;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=f0, g=g);
            basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes, normalization=Normalization.uMax);
            [FExpected, GExpected, hExpected] = testCase.v1ModesAtWavenumber(N0, b, zDomain, z, nModes, g, k);

            testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-11)
            testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
            testCase.verifyModeShapes(basisSet.G(z), GExpected, 1e-10)
            testCase.verifyModeShapes(basisSet.F(z), FExpected, 1e-10)
            testCase.verifyGreaterThan(basisSet.F(zDomain(2)), zeros(1,nModes))
        end

        function fixedFrequencyGModesMatchV1Oracle(testCase)
            [N0, b, zDomain, z, nModes, ~, g] = testCase.profile();
            omega = 0.8*N0;
            evp = IMEigenvalueProblem.waveModesAtFrequency(omega=omega, f0=0, g=g);
            basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes, normalization=Normalization.uMax);
            [FExpected, GExpected, hExpected] = testCase.v1ModesAtFrequency(N0, b, zDomain, z, nModes, g, omega);

            testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-11)
            testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
            testCase.verifyModeShapes(basisSet.G(z), GExpected, 1e-10)
            testCase.verifyModeShapes(basisSet.F(z), FExpected, 1e-10)
            testCase.verifyGreaterThan(basisSet.F(zDomain(2)), zeros(1,nModes))
        end

        function hydrostaticGModesMatchV1Oracle(testCase)
            [N0, b, zDomain, z, nModes, f0, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g);
            basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes, normalization=Normalization.uMax);
            [FExpected, GExpected, hExpected] = testCase.v1ModesAtFrequency(N0, b, zDomain, z, nModes, g, 0);

            testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-11)
            testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
            testCase.verifyModeShapes(basisSet.G(z), GExpected, 1e-10)
            testCase.verifyModeShapes(basisSet.F(z), FExpected, 1e-10)
            testCase.verifyGreaterThan(basisSet.F(zDomain(2)), zeros(1,nModes))
            testCase.verifyLessThan(max(abs(basisSet.G(zDomain(1)))), 1e-11)
            testCase.verifyLessThan(max(abs(basisSet.G(zDomain(2)))), 1e-11)
        end

        function freeSurfaceWavenumberGModesMatchV1Oracle(testCase)
            [N0, b, zDomain, z, nModes, f0, g] = testCase.profile();
            k_star = sqrt((N0*N0 - f0*f0)/(g*diff(zDomain)));
            for k = [0.1*k_star k_star 10*k_star]
                evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=f0, g=g, ...
                    surfaceBoundary=IMBoundary.free());
                basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                    zDomain=zDomain, nModes=nModes, normalization=Normalization.uMax);
                [FExpected, GExpected, hExpected] = testCase.v1FreeSurfaceModesAtWavenumber( ...
                    N0, b, zDomain, z, nModes, g, k);
                [F, G] = testCase.surfaceNormalizedModes(basisSet, z);

                testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-11)
                testCase.verifyEqual(basisSet.modeNumber, [-1 1:(nModes-1)])
                testCase.verifyModeShapes(G, GExpected, 1e-10)
                testCase.verifyModeShapes(F, FExpected, 1e-10)
                testCase.verifyFreeSurfaceBoundaryResiduals(basisSet, zDomain, 1e-10)
            end
        end

        function freeSurfaceFrequencyGModesMatchV1Oracle(testCase)
            [N0, b, zDomain, z, nModes, f0, g] = testCase.profile();
            for omega = [0.1*N0 0]
                if omega == 0
                    evp = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g, surfaceBoundary=IMBoundary.free());
                else
                    evp = IMEigenvalueProblem.waveModesAtFrequency(omega=omega, f0=f0, g=g, ...
                        surfaceBoundary=IMBoundary.free());
                end
                basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                    zDomain=zDomain, nModes=nModes, normalization=Normalization.uMax);
                [FExpected, GExpected, hExpected] = testCase.v1FreeSurfaceModesAtFrequency( ...
                    N0, b, zDomain, z, nModes, g, omega);
                [F, G] = testCase.surfaceNormalizedModes(basisSet, z);

                testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-11)
                testCase.verifyEqual(basisSet.modeNumber, [-1 1:(nModes-1)])
                testCase.verifyModeShapes(G, GExpected, 1e-10)
                testCase.verifyModeShapes(F, FExpected, 1e-10)
                testCase.verifyFreeSurfaceBoundaryResiduals(basisSet, zDomain, 1e-10)
            end
        end

        function highFrequencyFreeSurfaceGModesMatchV1Oracle(testCase)
            [N0, b, zDomain, z, nModes, f0, g] = testCase.profile();
            for omega = [N0 10*N0]
                evp = IMEigenvalueProblem.waveModesAtFrequency(omega=omega, f0=f0, g=g, ...
                    surfaceBoundary=IMBoundary.free());
                basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                    zDomain=zDomain, nModes=nModes, normalization=Normalization.surfacePressure);
                [FExpected, GExpected, hExpected] = testCase.v1FreeSurfaceModesAtFrequency( ...
                    N0, b, zDomain, z, nModes, g, omega);

                testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-11)
                testCase.verifyEqual(basisSet.modeNumber, -1)
                testCase.verifyModeShapes(basisSet.G(z), GExpected, 1e-10)
                testCase.verifyModeShapes(basisSet.F(z), FExpected, 1e-10)
                testCase.verifyFreeSurfaceBoundaryResiduals(basisSet, zDomain, 1e-10)
            end
        end

        function hydrostaticFModesReuseHydrostaticGBaroclinicBranches(testCase)
            [N0, b, zDomain, z, nModes, f0, g] = testCase.profile();
            evpF = IMEigenvalueProblem.hydrostaticFModes(g=g);
            evpG = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g);
            basisF = IMBasisSet.exponentialStratification(evp=evpF, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes, normalization=Normalization.geostrophic);
            basisG = IMBasisSet.exponentialStratification(evp=evpG, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes-1, normalization=Normalization.geostrophic);

            testCase.verifyEqual(basisF.modeNumber, 0:(nModes-1))
            testCase.verifyEqual(basisF.eigenvalues(1), 0, AbsTol=0)
            testCase.verifyTrue(isinf(basisF.h(1)))
            testCase.verifyEquivalentDepths(basisF.h(2:end), basisG.h, 1e-11)
            testCase.verifyGreaterThan(basisF.F(zDomain(2)), zeros(1,nModes))
            testCase.verifyLessThan(max(abs(basisF.G(zDomain(2)))), 1e-10)
            testCase.verifyLessThan(max(abs(basisF.G(zDomain(1)))), 1e-10)
            GF = basisF.G(z);
            FF = basisF.F(z);
            testCase.verifyLessThan(max(abs(GF(:,1))), 1e-12)
            testCase.verifyLessThan(max(abs(FF(:,1) - 1)), 1e-12)
            testCase.verifyModeShapes(GF(:,2:end), basisG.G(z), 1e-10)
            testCase.verifyModeShapes(FF(:,2:end), basisG.F(z), 1e-10)

            basisF.normalization = Normalization.unity;
            gramF = basisF.gramMatrix("F");
            testCase.verifyLessThan(max(abs(diag(gramF).' - 1)), 1e-8)

            basisF.normalization = Normalization.surfacePressure;
            testCase.verifyEqual(basisF.F(zDomain(2)), ones(1,nModes), AbsTol=1e-10)
        end

        function freeSurfaceSingleModeReturnsSurfaceBranch(testCase)
            [N0, b, zDomain, ~, ~, f0, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g, surfaceBoundary=IMBoundary.free());
            basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=1);

            testCase.verifyEqual(basisSet.modeNumber, -1)
            testCase.verifyGreaterThan(basisSet.h, diff(zDomain)/2)
            testCase.verifyFreeSurfaceBoundaryResiduals(basisSet, zDomain, 1e-10)
        end

        function surfacePressureNormalizationScalesExponentialGModes(testCase)
            [N0, b, zDomain, ~, nModes, f0, g] = testCase.profile();
            cases = {
                IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g), false
                IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g, ...
                    surfaceBoundary=IMBoundary.free()), true
            };

            for iCase = 1:size(cases,1)
                evp = cases{iCase,1};
                hasFreeSurface = cases{iCase,2};
                basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                    zDomain=zDomain, nModes=nModes, normalization=Normalization.surfacePressure);

                testCase.verifyEqual(basisSet.F(zDomain(2)), ones(1,nModes), AbsTol=1e-10)
                if hasFreeSurface
                    testCase.verifyEqual(basisSet.G(zDomain(2)), ones(1,nModes), AbsTol=1e-10)
                else
                    testCase.verifyLessThan(max(abs(basisSet.G(zDomain(2)))), 1e-10)
                end
            end
        end

        function allSolversMatchExponentialWavenumberAnalyticalBasis(testCase)
            [N0, b, zDomain, z, nModes, f0, g] = testCase.profile();
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            exactBasis = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes, normalization=Normalization.uMax);
            cases = testCase.solverCases(N2, zDomain);

            for iCase = 1:size(cases,1)
                solver = cases{iCase,2};
                hTolerance = cases{iCase,3};
                shapeTolerance = cases{iCase,4};
                gramTolerance = cases{iCase,5};
                basisSet = solver.solveEVP(evp, nModes=nModes);
                basisSet.normalization = Normalization.uMax;

                testCase.verifyEquivalentDepths(basisSet.h, exactBasis.h, hTolerance)
                testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
                testCase.verifyLessThan(testCase.relativeError(basisSet.G(z), exactBasis.G(z)), shapeTolerance)
                testCase.verifyLessThan(testCase.relativeError(basisSet.F(z), exactBasis.F(z)), shapeTolerance)
                testCase.verifyGreaterThan(basisSet.F(zDomain(2)), zeros(1,nModes))
                testCase.verifyLessThan(max(abs(basisSet.G(zDomain(1)))), shapeTolerance)
                testCase.verifyLessThan(max(abs(basisSet.G(zDomain(2)))), shapeTolerance)

                basisSet.normalization = Normalization.unity;
                gram = basisSet.gramMatrix("G");
                testCase.verifyLessThan(max(abs(diag(gram).' - 1)), gramTolerance)
                testCase.verifyLessThan(testCase.offDiagonalNorm(gram), gramTolerance)
            end
        end

        function allSolversMatchFreeSurfaceExponentialWavenumberAnalyticalBasis(testCase)
            [N0, b, zDomain, z, nModes, f0, g] = testCase.profile();
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g, ...
                surfaceBoundary=IMBoundary.free());
            exactBasis = IMBasisSet.exponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes, normalization=Normalization.uMax);
            [FExact, GExact] = testCase.surfaceNormalizedModes(exactBasis, z);
            cases = testCase.solverCases(N2, zDomain);

            for iCase = 1:size(cases,1)
                solver = cases{iCase,2};
                hTolerance = max(cases{iCase,3}, 3e-3);
                shapeTolerance = max(cases{iCase,4}, 3e-3);
                gramTolerance = max(cases{iCase,5}, 1e-4);
                basisSet = solver.solveEVP(evp, nModes=nModes);
                basisSet.normalization = Normalization.uMax;
                [F, G] = testCase.surfaceNormalizedModes(basisSet, z);

                testCase.verifyEquivalentDepths(basisSet.h, exactBasis.h, hTolerance)
                testCase.verifyEqual(basisSet.modeNumber, [-1 1:(nModes-1)])
                testCase.verifyLessThan(testCase.relativeError(G, GExact), shapeTolerance)
                testCase.verifyLessThan(testCase.relativeError(F, FExact), shapeTolerance)
                testCase.verifyFreeSurfaceBoundaryResiduals(basisSet, zDomain, shapeTolerance)

                basisSet.normalization = Normalization.unity;
                gram = basisSet.gramMatrix("G");
                testCase.verifyLessThan(max(abs(diag(gram).' - 1)), gramTolerance)
                testCase.verifyLessThan(testCase.offDiagonalNorm(gram), gramTolerance)
            end
        end

        function unsupportedAnalyticalBoundariesThrow(testCase)
            [N0, b, zDomain, ~, nModes, ~, g] = testCase.profile();
            evpFreeSurfaceF = IMEigenvalueProblem.hydrostaticFModes(g=g, surfaceBoundary=IMBoundary.free());
            evpRigidTooFast = IMEigenvalueProblem.waveModesAtFrequency(omega=N0, g=g);
            evpNoSlipBottom = IMEigenvalueProblem.hydrostaticGModes(g=g, bottomBoundary=IMBoundary.noSlip());
            evpFreeBottom = IMEigenvalueProblem.hydrostaticGModes(g=g, bottomBoundary=IMBoundary.free());
            customSurface = IMBoundary.custom(left=IMOperator().plus(derivativeOrder=0));
            evpCustomSurface = IMEigenvalueProblem.hydrostaticGModes(g=g, surfaceBoundary=customSurface);

            testCase.verifyError(@() IMBasisSet.exponentialStratification(evp=evpFreeSurfaceF, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetExponentialStratification:UnsupportedBoundary")
            testCase.verifyError(@() IMBasisSet.exponentialStratification(evp=evpRigidTooFast, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetExponentialStratification:UnsupportedFrequency")
            testCase.verifyError(@() IMBasisSet.exponentialStratification(evp=evpNoSlipBottom, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetExponentialStratification:UnsupportedBoundary")
            testCase.verifyError(@() IMBasisSet.exponentialStratification(evp=evpFreeBottom, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetExponentialStratification:UnsupportedBoundary")
            testCase.verifyError(@() IMBasisSet.exponentialStratification(evp=evpCustomSurface, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetExponentialStratification:UnsupportedBoundary")
        end

        function v2AnalyticalCodeAvoidsDeferredDependencyAndOldBoundaryTokens(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            files = [
                fullfile(repoRoot, "@IMBasisSetExponentialStratification", "IMBasisSetExponentialStratification.m")
                fullfile(repoRoot, "@IMBasisSetConstantStratification", "IMBasisSetConstantStratification.m")
            ];
            forbidden = ["cheb" + "fun", "rigid" + "Lid", "free" + "Slip"];
            for iFile = 1:length(files)
                source = string(fileread(files(iFile)));
                for iToken = 1:length(forbidden)
                    testCase.verifyFalse(contains(source, forbidden(iToken), IgnoreCase=true))
                end
            end
        end
    end

    methods (Access = private)
        function [N0, b, zDomain, z, nModes, f0, g] = profile(~)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            z = linspace(zDomain(1), zDomain(2), 96).';
            nModes = 4;
            f0 = 0;
            g = 9.81;
        end

        function [F, G, h] = v1ModesAtWavenumber(testCase, N0, b, zDomain, z, nModes, g, k)
            direct = testCase.v1Oracle(N0, b, zDomain, z, nModes, g);
            [F, G, h] = direct.modesAtWavenumber(k);
        end

        function [F, G, h] = v1ModesAtFrequency(testCase, N0, b, zDomain, z, nModes, g, omega)
            direct = testCase.v1Oracle(N0, b, zDomain, z, nModes, g);
            [F, G, h] = direct.modesAtFrequency(omega);
        end

        function [F, G, h] = v1FreeSurfaceModesAtWavenumber(testCase, N0, b, zDomain, z, nModes, g, k)
            direct = testCase.v1FreeSurfaceOracle(N0, b, zDomain, z, nModes, g);
            [F, G, h] = direct.modesAtWavenumber(k);
        end

        function [F, G, h] = v1FreeSurfaceModesAtFrequency(testCase, N0, b, zDomain, z, nModes, g, omega)
            direct = testCase.v1FreeSurfaceOracle(N0, b, zDomain, z, nModes, g);
            [F, G, h] = direct.modesAtFrequency(omega);
        end

        function direct = v1Oracle(~, N0, b, zDomain, z, nModes, g)
            direct = InternalModesExponentialStratification(N0=N0, b=b, zIn=zDomain, ...
                zOut=z, latitude=0, nModes=nModes, g=g);
            direct.normalization = Normalization.uMax;
        end

        function direct = v1FreeSurfaceOracle(testCase, N0, b, zDomain, z, nModes, g)
            direct = testCase.v1Oracle(N0, b, zDomain, z, nModes, g);
            direct.upperBoundary = UpperBoundary.freeSurface;
        end

        function cases = solverCases(testCase, N2, zDomain)
            nSpectral = 160;
            nFiniteDifference = 240;
            cases = {
                "spectral", IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nSpectral), 1e-8, 1e-7, 1e-8
                "wkb", IMSolverWKBSpectral(N2=N2, zDomain=zDomain, nEVP=nSpectral), 1e-5, 1e-5, 1e-5
                "density", IMSolverDensitySpectral(N2=N2, zDomain=zDomain, nEVP=nSpectral), 1e-4, 2e-3, 1e-4
                "finiteDifference", testCase.finiteDifferenceSolver(N2, zDomain, nFiniteDifference), 1e-5, 1e-3, 1e-4
            };
        end

        function solver = finiteDifferenceSolver(~, N2, zDomain, nGrid)
            z = linspace(zDomain(1), zDomain(2), nGrid).';
            solver = IMSolverFiniteDifference(z=z, N2=N2);
        end

        function verifyModeShapes(testCase, actual, expected, tolerance)
            testCase.verifyLessThan(testCase.relativeError(actual, expected), tolerance)
        end

        function verifyEquivalentDepths(testCase, actual, expected, relativeTolerance)
            actual = reshape(actual,1,[]);
            expected = reshape(expected,1,[]);
            testCase.verifyEqual(actual, expected, RelTol=relativeTolerance, AbsTol=relativeTolerance)
        end

        function value = relativeError(~, actual, expected)
            value = norm(actual - expected, "fro")/max(1,norm(expected, "fro"));
        end

        function [F, G] = surfaceNormalizedModes(~, basisSet, z)
            FSurface = basisSet.F(basisSet.zDomain(2));
            F = basisSet.F(z)./FSurface;
            G = basisSet.G(z)./FSurface;
        end

        function verifyFreeSurfaceBoundaryResiduals(testCase, basisSet, zDomain, tolerance)
            testCase.verifyGreaterThan(basisSet.F(zDomain(2)), zeros(size(basisSet.F(zDomain(2)))))
            testCase.verifyLessThan(max(abs(basisSet.F(zDomain(2)) - basisSet.G(zDomain(2)))), tolerance)
            testCase.verifyLessThan(max(abs(basisSet.G(zDomain(1)))), tolerance)
        end

        function value = offDiagonalNorm(~, matrix)
            value = norm(matrix - diag(diag(matrix)), "fro");
        end
    end
end
