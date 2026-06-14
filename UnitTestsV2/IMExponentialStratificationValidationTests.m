classdef IMExponentialStratificationValidationTests < matlab.unittest.TestCase

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
        function rigidHydrostaticGModesEvaluateExponentialBasis(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);

            basisSet = IMBasisSetExponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 24).';

            testCase.verifyClass(basisSet, "IMBasisSetExponentialStratification")
            testCase.verifyTrue(all(isfinite(basisSet.h)))
            testCase.verifySize(basisSet.G(z), [24 3])
            testCase.verifySize(basisSet.F(z), [24 3])
            testCase.verifyEqual(basisSet.modeNumber, 1:3)
        end

        function freeSurfaceHydrostaticGModesIncludeBoundaryMode(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            freeSurface = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, surfaceBoundary=freeSurface);

            basisSet = IMBasisSetExponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=3);

            testCase.verifyEqual(basisSet.modeNumber(1), -1)
            testCase.verifyTrue(all(isfinite(basisSet.roots)))
        end

        function hydrostaticFModesIncludeNullBranch(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain);

            basisSet = IMBasisSetExponentialStratification(evp=evp, N0=N0, b=b, ...
                zDomain=zDomain, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 12).';

            testCase.verifyEqual(basisSet.modeNumber(1), 0)
            F = basisSet.F(z, normalization=Normalization.unity);
            G = basisSet.G(z, normalization=Normalization.unity);
            factors = basisSet.normalizationFactors(Normalization.unity);
            testCase.verifyEqual(F(:,1), ones(size(z))/factors(1), AbsTol=1e-12)
            testCase.verifyEqual(G(:,1), zeros(size(z)), AbsTol=1e-12)
        end

        function analyticalFactoryCreatesExponentialBasis(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.waveModesAtFrequency(N2=N2, zDomain=zDomain, omega=1e-3);

            basisSet = IMBasisSet.exponentialStratification(evp=evp, N0=N0, ...
                b=b, zDomain=zDomain, nModes=2);

            testCase.verifyClass(basisSet, "IMBasisSetExponentialStratification")
            testCase.verifyEqual(basisSet.evp.parameters.omega, 1e-3, AbsTol=0)
        end

        function providedEVPDomainMustMatchAnalyticalDomain(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);

            testCase.verifyError(@() IMBasisSetExponentialStratification(evp=evp, ...
                N0=N0, b=b, zDomain=[-4000 0], nModes=2), ...
                "IMBasisSetExponentialStratification:DomainMismatch")
        end

        function unsupportedBoundaryIsRejected(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-5000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, bottomBoundary=IMBoundaryCondition.neumann());

            testCase.verifyError(@() IMBasisSetExponentialStratification(evp=evp, ...
                N0=N0, b=b, zDomain=zDomain, nModes=2), ...
                "IMBasisSetExponentialStratification:UnsupportedBoundary")
        end
    end
end
