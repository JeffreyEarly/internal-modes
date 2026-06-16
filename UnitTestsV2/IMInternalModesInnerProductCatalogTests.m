classdef IMInternalModesInnerProductCatalogTests < matlab.unittest.TestCase

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
        function simpleHydrostaticDiagnosticRowsAreKnown(testCase)
            [N2, zDomain] = testCase.profile();
            gEVP = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
            fEVP = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain);

            fSpec = gEVP.innerProduct("F");
            gSpec = fEVP.innerProduct("G");

            testCase.verifyTrue(fSpec.hasInnerProduct)
            testCase.verifyTrue(gSpec.hasInnerProduct)
            testCase.verifyEmpty(fSpec.endpointInnerProductTerms)
            testCase.verifyEmpty(gSpec.endpointInnerProductTerms)
            testCase.verifyFalse(isfield(fSpec, "status"))
            testCase.verifyFalse(isfield(fSpec, "endpointFunctionals"))
            testCase.verifyTrue(contains(fSpec.reason, "interior-only"))
        end

        function diagnosticEndpointValueRowsAreKnown(testCase)
            [N2, zDomain, ~, g] = testCase.profile();
            gEVP = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=3));
            fRobinEVP = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, g=g, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=3));
            fLinkedEVP = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, g=g, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=4, c=1, d=0));

            fSpec = gEVP.innerProduct("F");
            gRobinSpec = fRobinEVP.innerProduct("G");
            gLinkedSpec = fLinkedEVP.innerProduct("G");

            testCase.verifyTrue(fSpec.hasInnerProduct)
            testCase.verifyEqual(fSpec.endpointInnerProductTerms(1).catalogCase, "G-P5")
            testCase.verifyEqual(fSpec.endpointInnerProductTerms(1).variable, "F")
            testCase.verifyEqual(fSpec.endpointInnerProductTerms(1).coefficient, -3/2, AbsTol=0)
            testCase.verifyTrue(gRobinSpec.hasInnerProduct)
            testCase.verifyEqual(gRobinSpec.endpointInnerProductTerms(1).catalogCase, "F-P2")
            testCase.verifyEqual(gRobinSpec.endpointInnerProductTerms(1).variable, "G")
            testCase.verifyEqual(gRobinSpec.endpointInnerProductTerms(1).coefficient, -3/(2*g), RelTol=1e-12)
            testCase.verifyTrue(gLinkedSpec.hasInnerProduct)
            testCase.verifyEqual(gLinkedSpec.endpointInnerProductTerms(1).catalogCase, "F-T4")
            testCase.verifyEqual(gLinkedSpec.endpointInnerProductTerms(1).variable, "F")
            testCase.verifyEqual(gLinkedSpec.endpointInnerProductTerms(1).coefficient, -g/2, RelTol=1e-12)
        end

        function endpointInnerProductTermsEnterGramMatrix(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=3));
            basisSet = solver.solveEVP(evp, nModes=2);
            spec = evp.innerProduct("F");

            z = solver.configuredForEVP(evp).innerProductGrid(zDomain);
            F = basisSet.F(z);
            gramInterior = zeros(2,2);
            for iMode = 1:2
                for jMode = iMode:2
                    value = basisSet.solver.integrateInnerProduct(z, F(:,iMode).*F(:,jMode), zDomain);
                    gramInterior(iMode,jMode) = value;
                    gramInterior(jMode,iMode) = value;
                end
            end
            FSurface = basisSet.F(zDomain(2));
            expected = gramInterior + spec.endpointInnerProductTerms(1).coefficient*(FSurface.'*FSurface);
            partialBounds = [zDomain(1) mean(zDomain)];
            zPartial = solver.configuredForEVP(evp).innerProductGrid(partialBounds);
            FPartial = basisSet.F(zPartial);
            expectedPartial = zeros(2,2);
            for iMode = 1:2
                for jMode = iMode:2
                    value = basisSet.solver.integrateInnerProduct(zPartial, FPartial(:,iMode).*FPartial(:,jMode), partialBounds);
                    expectedPartial(iMode,jMode) = value;
                    expectedPartial(jMode,iMode) = value;
                end
            end

            testCase.verifyEqual(basisSet.gramMatrix("F"), expected, RelTol=1e-10, AbsTol=1e-10)
            testCase.verifyEqual(basisSet.partialGramMatrix("F", partialBounds(1), partialBounds(2)), ...
                expectedPartial, RelTol=1e-10, AbsTol=1e-10)
        end

        function linkedEndpointInnerProductTermsEnterGramMatrix(testCase)
            [N2, zDomain, nEVP, g] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, g=g, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=4, c=1, d=0));
            basisSet = solver.solveEVP(evp, nModes=2);
            spec = evp.innerProduct("G");

            z = solver.configuredForEVP(evp).innerProductGrid(zDomain);
            G = basisSet.G(z);
            weight = N2(z)/g;
            gramInterior = zeros(2,2);
            for iMode = 1:2
                for jMode = iMode:2
                    value = basisSet.solver.integrateInnerProduct(z, weight.*G(:,iMode).*G(:,jMode), zDomain);
                    gramInterior(iMode,jMode) = value;
                    gramInterior(jMode,iMode) = value;
                end
            end
            FSurface = basisSet.F(zDomain(2));
            expected = gramInterior + spec.endpointInnerProductTerms(1).coefficient*(FSurface.'*FSurface);

            testCase.verifyEqual(basisSet.gramMatrix("G"), expected, RelTol=1e-10, AbsTol=1e-10)
        end

        function unsupportedHydrostaticRowsHaveNoInnerProduct(testCase)
            [N2, zDomain] = testCase.profile();
            gT3 = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=1, c=1, d=0));
            gT1 = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=0, c=0, d=1));
            fT2 = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=0, c=1, d=1));
            allFour = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=1, c=1, d=1));

            testCase.verifyUnavailable(gT3.innerProduct("F"), "derivative endpoint terms")
            testCase.verifyUnavailable(gT1.innerProduct("F"), "derivative endpoint terms")
            testCase.verifyUnavailable(fT2.innerProduct("G"), "derivative endpoint terms")
            testCase.verifyUnavailable(allFour.innerProduct("F"), "outside the value-only")
        end

        function unavailableDiagnosticGramMatrixThrows(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=0, c=0, d=1));
            basisSet = solver.solveEVP(evp, nModes=2);

            testCase.verifyError(@() basisSet.gramMatrix("F"), "IMInternalModesBasis:UnavailableInnerProduct")
        end

        function waveDiagnosticInnerProductsRemainUnavailable(testCase)
            [N2, zDomain, nEVP, g] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.waveModesAtFrequency(N2=N2, zDomain=zDomain, omega=1e-3, g=g);
            basisSet = solver.solveEVP(evp, nModes=2);
            spec = evp.innerProduct("F");

            testCase.verifyFalse(spec.hasInnerProduct)
            testCase.verifyTrue(contains(spec.reason, "geostrophic"))
            testCase.verifyError(@() basisSet.gramMatrix("F"), "IMInternalModesBasis:UnavailableInnerProduct")
        end

        function customGeostrophicFamilyUsesCatalogWithoutHydrostaticName(testCase)
            [N2, zDomain, ~, g] = testCase.profile();
            surface = IMBoundaryCondition(a=2, b=3);
            bottom = IMBoundaryCondition(a=-4, b=5);
            evp = IMInternalModes(name="unforced-APV-modes", formulation="F", ...
                modeFamily="geostrophic", N2=N2, zDomain=zDomain, ...
                p=@(z,ctx) 1./ctx.N2(z), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ones(size(z))/ctx.g, g=g, ...
                surfaceBoundary=surface, bottomBoundary=bottom);

            spec = evp.innerProduct("G");

            testCase.verifyEqual(evp.name, "unforced-APV-modes")
            testCase.verifyEqual(evp.modeFamily, "geostrophic")
            testCase.verifyTrue(spec.hasInnerProduct)
            testCase.verifyEqual([spec.endpointInnerProductTerms.catalogCase], ["F-P2" "F-P2"])
            testCase.verifyEqual(spec.endpointInnerProductTerms(1).coefficient, -3/(2*g), RelTol=1e-12)
            testCase.verifyEqual(spec.endpointInnerProductTerms(2).coefficient, -5/(4*g), RelTol=1e-12)
        end

        function geostrophicNormalizationCouplesFAndGForGFormulatedModes(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
            basisSet = solver.solveEVP(evp, nModes=2);
            basisSet.normalization = "geostrophic";

            testCase.verifyCoupledGeostrophicNormalization(basisSet)
        end

        function geostrophicNormalizationCouplesFAndGForFFormulatedModes(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain);
            basisSet = solver.solveEVP(evp, nModes=3);
            basisSet.normalization = "geostrophic";

            testCase.verifyCoupledGeostrophicNormalization(basisSet)
        end

        function geostrophicNormalizationUsesEndpointValueDiagnosticTerms(testCase)
            [N2, zDomain, nEVP, g] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, g=g, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=3));
            basisSet = solver.solveEVP(evp, nModes=2);
            basisSet = basisSet.addNormalization("raw", @(~,~) 1);
            basisSet.normalization = "raw";
            rawGGram = basisSet.gramMatrix("G");
            rawFactor = sqrt(abs(rawGGram(1,1)));
            factors = basisSet.normalizationFactors("geostrophic");

            testCase.verifyEqual( ...
                factors(1), ...
                rawFactor, RelTol=2e-8, AbsTol=2e-8)
        end
    end

    methods
        function verifyUnavailable(testCase, spec, reasonText)
            testCase.verifyFalse(spec.hasInnerProduct)
            testCase.verifyTrue(contains(spec.reason, reasonText))
            testCase.verifyFalse(isfield(spec, "status"))
        end

        function verifyCoupledGeostrophicNormalization(testCase, basisSet)
            baroclinic = basisSet.modeNumber ~= 0;
            gramG = basisSet.gramMatrix("G");
            gramF = basisSet.gramMatrix("F");

            testCase.verifyEqual(diag(gramG(baroclinic,baroclinic)).', ...
                ones(1,nnz(baroclinic)), RelTol=2e-8, AbsTol=2e-8)
            testCase.verifyEqual(diag(gramF(baroclinic,baroclinic)).', ...
                basisSet.h(baroclinic), RelTol=2e-7, AbsTol=2e-7)
        end
    end

    methods (Static, Access = private)
        function [N2, zDomain, nEVP, g] = profile()
            N0 = 5.2e-3;
            zDomain = [-1000 0];
            nEVP = 24;
            g = 9.81;
            N2 = @(z) N0*N0*(1 + 0.1*z/abs(zDomain(1)));
        end
    end
end
