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

        function unsupportedAnalyticalBoundariesThrow(testCase)
            [N0, b, zDomain, ~, nModes, ~, g] = testCase.profile();
            evpF = IMEigenvalueProblem.hydrostaticFModes(g=g);
            evpFreeSurface = IMEigenvalueProblem.hydrostaticGModes(g=g, surfaceBoundary=IMBoundary.free());
            evpNoSlipBottom = IMEigenvalueProblem.hydrostaticGModes(g=g, bottomBoundary=IMBoundary.noSlip());
            evpTooFast = IMEigenvalueProblem.waveModesAtFrequency(omega=N0, g=g);

            testCase.verifyError(@() IMBasisSet.exponentialStratification(evp=evpF, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetExponentialStratification:UnsupportedEVP")
            testCase.verifyError(@() IMBasisSet.exponentialStratification(evp=evpFreeSurface, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetExponentialStratification:UnsupportedBoundary")
            testCase.verifyError(@() IMBasisSet.exponentialStratification(evp=evpNoSlipBottom, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetExponentialStratification:UnsupportedBoundary")
            testCase.verifyError(@() IMBasisSet.exponentialStratification(evp=evpTooFast, N0=N0, b=b, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetExponentialStratification:UnsupportedFrequency")
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

        function direct = v1Oracle(~, N0, b, zDomain, z, nModes, g)
            direct = InternalModesExponentialStratification(N0=N0, b=b, zIn=zDomain, ...
                zOut=z, latitude=0, nModes=nModes, g=g);
            direct.normalization = Normalization.uMax;
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

        function value = offDiagonalNorm(~, matrix)
            value = norm(matrix - diag(diag(matrix)), "fro");
        end
    end
end
