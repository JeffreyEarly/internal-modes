classdef IMConstantStratificationValidationTests < matlab.unittest.TestCase

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
        function fixedWavenumberGModesMatchRigidBoundaryFormula(testCase)
            [N0, zDomain, z, nModes, f0, g] = testCase.profile();
            k = 1e-4;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=f0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, ...
                nModes=nModes, normalization=Normalization.wMax);
            [GExpected, FExpected, hExpected] = testCase.fixedWavenumberGFormula(z, zDomain, N0, g, nModes, k, f0);

            testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-12)
            testCase.verifyEqual(basisSet.eigenvalues, 1./hExpected, RelTol=1e-12)
            testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
            testCase.verifyEqual(basisSet.G(z), GExpected, AbsTol=1e-12)
            testCase.verifyEqual(basisSet.F(z), FExpected, AbsTol=1e-12)
        end

        function fixedFrequencyGModesMatchRigidBoundaryFormula(testCase)
            [N0, zDomain, z, nModes, ~, g] = testCase.profile();
            omega = 0.8*N0;
            evp = IMEigenvalueProblem.waveModesAtFrequency(omega=omega, f0=0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, ...
                nModes=nModes, normalization=Normalization.wMax);
            [GExpected, FExpected, hExpected] = testCase.fixedFrequencyGFormula(z, zDomain, N0, g, nModes, omega);

            testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-12)
            testCase.verifyEqual(basisSet.eigenvalues, 1./hExpected, RelTol=1e-12)
            testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
            testCase.verifyEqual(basisSet.G(z), GExpected, AbsTol=1e-12)
            testCase.verifyEqual(basisSet.F(z), FExpected, AbsTol=1e-12)
        end

        function hydrostaticGModesMatchRigidBoundaryFormula(testCase)
            [N0, zDomain, z, nModes, f0, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, ...
                nModes=nModes, normalization=Normalization.wMax);
            [GExpected, FExpected, hExpected] = testCase.fixedFrequencyGFormula(z, zDomain, N0, g, nModes, 0);

            testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-12)
            testCase.verifyEqual(basisSet.eigenvalues, 1./hExpected, RelTol=1e-12)
            testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
            testCase.verifyEqual(basisSet.G(z), GExpected, AbsTol=1e-12)
            testCase.verifyEqual(basisSet.F(z), FExpected, AbsTol=1e-12)
        end

        function hydrostaticFModesIncludeExactNullMode(testCase)
            [N0, zDomain, z, nModes, ~, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, ...
                nModes=nModes, normalization=Normalization.uMax);
            [GExpected, FExpected, dFdzExpected, hExpected, eigenvaluesExpected] = ...
                testCase.hydrostaticFFormula(z, zDomain, N0, g, nModes);

            testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-12)
            testCase.verifyEqual(basisSet.eigenvalues, eigenvaluesExpected, AbsTol=1e-12)
            testCase.verifyEqual(basisSet.modeNumber, 0:(nModes-1))
            testCase.verifyEqual(basisSet.G(z), GExpected, AbsTol=1e-12)
            testCase.verifyEqual(basisSet.F(z), FExpected, AbsTol=1e-12)
            testCase.verifyEqual(GExpected, -(g/(N0*N0))*dFdzExpected, AbsTol=1e-12)
        end

        function hydrostaticFGeostrophicGramMatchesExactMetrics(testCase)
            [N0, zDomain, ~, nModes, ~, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);
            basisSet.normalization = Normalization.geostrophic;

            gramG = basisSet.gramMatrix("G");
            gramF = basisSet.gramMatrix("F");

            testCase.verifyEqual(gramG(1,1), 0, AbsTol=1e-12)
            testCase.verifyEqual(gramF(1,1), diff(zDomain), AbsTol=1e-12)
            testCase.verifyEqual(diag(gramG(2:end,2:end)).', ones(1,nModes-1), AbsTol=1e-12)
            testCase.verifyEqual(diag(gramF(2:end,2:end)).', basisSet.h(2:end), RelTol=1e-12)
            testCase.verifyLessThan(testCase.offDiagonalNorm(gramG), 1e-12)
            testCase.verifyLessThan(testCase.offDiagonalNorm(gramF), 1e-12)
        end

        function freeSurfaceGModesMatchV1OracleAndBoundaryResidual(testCase)
            [N0, zDomain, z, nModes, ~, g] = testCase.profile();
            k = 1e-4;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=0, g=g, ...
                surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.rigid());
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);

            direct = InternalModesConstantStratification(N0=N0, zIn=zDomain, zOut=z, latitude=0, nModes=nModes, g=g);
            direct.upperBoundary = UpperBoundary.freeSurface;
            direct.normalization = Normalization.kConstant;
            [FExpected, GExpected, hExpected] = direct.modesAtWavenumber(k);

            testCase.verifyEquivalentDepths(basisSet.h, hExpected, 1e-10)
            testCase.verifyEqual(basisSet.modeNumber, [-1 1:(nModes-1)])
            testCase.verifyEqual(basisSet.G(z), GExpected, AbsTol=1e-8)
            testCase.verifyEqual(basisSet.F(z), FExpected, AbsTol=1e-8)
            testCase.verifyLessThan(max(abs(basisSet.G(zDomain(1)))), 1e-10)
            testCase.verifyLessThan(max(abs(basisSet.F(zDomain(2)) - basisSet.G(zDomain(2)))), 1e-8)
        end

        function unsupportedHydrostaticFAnalyticalBoundariesThrow(testCase)
            [N0, zDomain, ~, nModes, ~, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g, surfaceBoundary=IMBoundary.noSlip(), ...
                bottomBoundary=IMBoundary.rigid());

            testCase.verifyError(@() IMBasisSet.constantStratification(evp=evp, N0=N0, ...
                zDomain=zDomain, nModes=nModes), "IMBasisSetConstantStratification:UnsupportedBoundary")
        end

        function allSolversMatchConstantWavenumberAnalyticalBasis(testCase)
            [N0, zDomain, z, nModes, f0, g] = testCase.profile();
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            exactBasis = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, ...
                nModes=nModes, normalization=Normalization.wMax);
            cases = testCase.solverCases(N2, zDomain);

            for iCase = 1:size(cases,1)
                solver = cases{iCase,2};
                hTolerance = cases{iCase,3};
                shapeTolerance = cases{iCase,4};
                gramTolerance = cases{iCase,5};
                basisSet = solver.solveEVP(evp, nModes=nModes);
                basisSet.normalization = Normalization.wMax;
                G = basisSet.G(z);
                F = basisSet.F(z);
                FSurface = basisSet.F(zDomain(2));

                testCase.verifyEquivalentDepths(basisSet.h, exactBasis.h, hTolerance)
                testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
                testCase.verifyLessThan(norm(G - exactBasis.G(z), "fro")/norm(exactBasis.G(z), "fro"), shapeTolerance)
                testCase.verifyLessThan(norm(F - exactBasis.F(z), "fro")/norm(exactBasis.F(z), "fro"), shapeTolerance)
                testCase.verifyTrue(all(FSurface > 0))
                testCase.verifyLessThan(max(abs(basisSet.G(zDomain(1)))), shapeTolerance)
                testCase.verifyLessThan(max(abs(basisSet.G(zDomain(2)))), shapeTolerance)

                basisSet.normalization = Normalization.unity;
                gram = basisSet.gramMatrix("G");
                testCase.verifyLessThan(max(abs(diag(gram).' - 1)), gramTolerance)
                testCase.verifyLessThan(testCase.offDiagonalNorm(gram), gramTolerance)
            end
        end

        function allSolversMatchConstantHydrostaticFAnalyticalBasis(testCase)
            [N0, zDomain, z, nModes, ~, g] = testCase.profile();
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);
            exactBasis = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, ...
                nModes=nModes, normalization=Normalization.uMax);
            cases = testCase.solverCases(N2, zDomain);

            for iCase = 1:size(cases,1)
                solver = cases{iCase,2};
                hTolerance = cases{iCase,3};
                shapeTolerance = cases{iCase,4};
                basisSet = solver.solveEVP(evp, nModes=nModes);
                basisSet.normalization = Normalization.uMax;
                G = basisSet.G(z);
                F = basisSet.F(z);
                FSurface = basisSet.F(zDomain(2));

                testCase.verifyLessThan(abs(basisSet.eigenvalues(1)), 1e-8)
                testCase.verifyEquivalentDepths(basisSet.h(2:end), exactBasis.h(2:end), hTolerance)
                testCase.verifyEqual(basisSet.modeNumber, 0:(nModes-1))
                GError = norm(G - exactBasis.G(z), "fro")/max(1,norm(exactBasis.G(z), "fro"));
                testCase.verifyLessThan(GError, shapeTolerance)
                testCase.verifyLessThan(norm(F - exactBasis.F(z), "fro")/norm(exactBasis.F(z), "fro"), shapeTolerance)
                testCase.verifyTrue(all(FSurface > 0))
                testCase.verifyLessThan(max(abs(G(:,1))), shapeTolerance)
                testCase.verifyLessThan(max(abs(F(:,1) - 1)), shapeTolerance)
            end
        end
    end

    methods (Access = private)
        function [N0, zDomain, z, nModes, f0, g] = profile(~)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            z = linspace(zDomain(1), zDomain(2), 96).';
            nModes = 5;
            f0 = 1e-4;
            g = 9.81;
        end

        function [G, F, h] = fixedWavenumberGFormula(~, z, zDomain, N0, g, nModes, k, f0)
            D = diff(zDomain);
            modeNumber = 1:nModes;
            m = modeNumber*pi/D;
            h = (N0*N0 - f0*f0)./(g*(k*k + m.*m));
            [G, F] = IMConstantStratificationValidationTests.rigidGVariables(z, zDomain, modeNumber, m, h);
        end

        function [G, F, h] = fixedFrequencyGFormula(~, z, zDomain, N0, g, nModes, omega)
            D = diff(zDomain);
            modeNumber = 1:nModes;
            m = modeNumber*pi/D;
            h = (N0*N0 - omega*omega)./(g*m.*m);
            [G, F] = IMConstantStratificationValidationTests.rigidGVariables(z, zDomain, modeNumber, m, h);
        end

        function [G, F, dFdz, h, eigenvalues] = hydrostaticFFormula(~, z, zDomain, N0, g, nModes)
            D = diff(zDomain);
            modeNumber = 0:(nModes-1);
            m = modeNumber*pi/D;
            h = [Inf (N0*N0)./(g*m(2:end).*m(2:end))];
            eigenvalues = [0 1./h(2:end)];
            s = z(:) - zDomain(1);
            signValue = (-1).^modeNumber;
            F = cos(s*m).*signValue;
            dFdz = -sin(s*m).*(signValue.*m);
            G = -(g/(N0*N0))*dFdz;
            G(:,1) = 0;
        end

        function cases = solverCases(testCase, N2, zDomain)
            nSpectral = 64;
            nFiniteDifference = 160;
            cases = {
                "spectral", IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nSpectral), 1e-8, 1e-3, 1e-10
                "wkb", IMSolverWKBSpectral(N2=N2, zDomain=zDomain, nEVP=nSpectral), 1e-8, 1e-3, 1e-10
                "density", IMSolverDensitySpectral(N2=N2, zDomain=zDomain, nEVP=nSpectral), 1e-8, 1e-3, 1e-10
                "finiteDifference", testCase.finiteDifferenceSolver(N2, zDomain, nFiniteDifference), 1e-5, 1e-3, 1e-6
            };
        end

        function verifyEquivalentDepths(testCase, actual, expected, relativeTolerance)
            actual = reshape(actual,1,[]);
            expected = reshape(expected,1,[]);
            finiteValues = isfinite(actual) & isfinite(expected);
            testCase.verifyEqual(actual(finiteValues), expected(finiteValues), ...
                RelTol=relativeTolerance, AbsTol=relativeTolerance)
            testCase.verifyEqual(isinf(abs(actual)), isinf(abs(expected)))
        end

        function value = offDiagonalNorm(~, matrix)
            value = norm(matrix - diag(diag(matrix)), "fro");
        end

        function solver = finiteDifferenceSolver(~, N2, zDomain, nGrid)
            z = linspace(zDomain(1), zDomain(2), nGrid).';
            solver = IMSolverFiniteDifference(z=z, N2=N2);
        end
    end

    methods (Static, Access = private)
        function [G, F] = rigidGVariables(z, zDomain, modeNumber, m, h)
            s = z(:) - zDomain(1);
            signValue = (-1).^modeNumber;
            G = sin(s*m).*signValue;
            F = cos(s*m).*(signValue.*h.*m);
        end
    end
end
