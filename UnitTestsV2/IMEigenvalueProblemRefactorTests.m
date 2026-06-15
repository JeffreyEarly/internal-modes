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
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMEigenvalueProblem(name="constant", zDomain=zDomain, p=2, q=3, r=4);
            configuredSolver = solver.configuredForEVP(evp);

            [A, B] = evp.assemble(solver);
            expectedA = -2*configuredSolver.physicalDerivativeMatrix(2) + 3*configuredSolver.physicalDerivativeMatrix(0);
            expectedB = 4*configuredSolver.physicalDerivativeMatrix(0);
            interiorRows = 2:(nEVP-1);

            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-11)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-11)
            testCase.verifyEqual(A(1,:), -configuredSolver.T(1,:), AbsTol=1e-11)
            testCase.verifyEqual(B(1,:), zeros(1,nEVP), AbsTol=1e-11)
            testCase.verifyEqual(A(end,:), -configuredSolver.T(end,:), AbsTol=1e-11)
            testCase.verifyEqual(B(end,:), zeros(1,nEVP), AbsTol=1e-11)
        end

        function canonicalAssemblyIncludesGridDerivativeOfVariableP(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            p = @(z,~) 1 + (z/1000).^2;
            evp = IMEigenvalueProblem(name="variableP", zDomain=zDomain, p=p, q=0, r=1);
            configuredSolver = solver.configuredForEVP(evp);

            [A, B] = evp.assemble(solver);
            pValues = p(configuredSolver.zNative, struct());
            pzValues = configuredSolver.differentiateGridValues(pValues, 1);
            expectedA = -diag(pValues)*configuredSolver.physicalDerivativeMatrix(2) ...
                - diag(pzValues)*configuredSolver.physicalDerivativeMatrix(1);
            expectedB = configuredSolver.physicalDerivativeMatrix(0);
            interiorRows = 2:(nEVP-1);

            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-9)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-11)
        end

        function endpointRowsUseCanonicalBoundaryCoefficients(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            surface = IMBoundaryCondition(a=2, b=3);
            bottom = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMEigenvalueProblem(zDomain=zDomain, p=5, surfaceBoundary=surface, bottomBoundary=bottom);
            configuredSolver = solver.configuredForEVP(evp);

            [A, B] = evp.assemble(solver);
            D0 = configuredSolver.physicalDerivativeMatrix(0);
            D1 = configuredSolver.physicalDerivativeMatrix(1);

            testCase.verifyEqual(A(1,:), -2*D0(1,:) + 15*D1(1,:), AbsTol=1e-11)
            testCase.verifyEqual(B(1,:), zeros(1,nEVP), AbsTol=1e-11)
            testCase.verifyEqual(A(end,:), 5*D1(end,:), AbsTol=1e-11)
            testCase.verifyEqual(B(end,:), D0(end,:), AbsTol=1e-11)
        end

        function internalModeFactoriesAssembleCanonicalForms(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            k = 1e-4;
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.waveModesAtWavenumber(N2=N2, zDomain=zDomain, k=k, f0=f0, g=g);
            configuredSolver = solver.configuredForEVP(evp);

            [A, B] = evp.assemble(solver);
            expectedA = -configuredSolver.physicalDerivativeMatrix(2) + k*k*configuredSolver.physicalDerivativeMatrix(0);
            expectedB = diag((N2(configuredSolver.zNative) - f0*f0)/g)*configuredSolver.physicalDerivativeMatrix(0);
            interiorRows = 2:(nEVP-1);

            testCase.verifyClass(evp, "IMInternalModes")
            testCase.verifyEqual(evp.formulation, "G")
            testCase.verifyEqual(evp.parameters.k, k, AbsTol=0)
            testCase.verifyEqual(evp.parameters.f0, f0, AbsTol=0)
            testCase.verifyEqual(evp.parameters.g, g, AbsTol=0)
            testCase.verifyEqual(evp.parameters.formulation, "G")
            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-11)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-11)
        end

        function evpParametersEnterCoefficientContext(testCase)
            [~, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMEigenvalueProblem(name="parameterized", zDomain=zDomain, ...
                p=@(z,ctx) ctx.alpha*ones(size(z)), q=0, r=1, ...
                parameters=struct("alpha",2));
            configuredSolver = solver.configuredForEVP(evp);

            context = evp.contextForSolver(configuredSolver);
            [A, ~] = evp.assemble(solver);
            expectedA = -2*configuredSolver.physicalDerivativeMatrix(2);
            interiorRows = 2:(nEVP-1);

            testCase.verifyEqual(evp.parameters.alpha, 2, AbsTol=0)
            testCase.verifyEqual(context.alpha, 2, AbsTol=0)
            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-11)
        end

        function innerProductSpecsUseCanonicalAndPhysicalVariables(testCase)
            [N2, zDomain, ~, ~, g] = testCase.profile();
            canonicalEVP = IMEigenvalueProblem(zDomain=zDomain, p=2, q=0, r=3);
            canonicalSpec = canonicalEVP.innerProduct();
            internalEVP = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, g=g);
            defaultSpec = internalEVP.innerProduct();
            gSpec = internalEVP.innerProduct("G");
            fSpec = internalEVP.innerProduct("F");
            z = linspace(zDomain(1), zDomain(2), 5).';
            context.N2 = N2;
            context.g = g;
            removedBoundaryFlag = "has" + "KnownBoundaryWeights";

            testCase.verifyEqual(canonicalSpec.variable, "u")
            testCase.verifyEqual(canonicalSpec.interiorWeight, canonicalEVP.r)
            testCase.verifyFalse(isfield(canonicalSpec, removedBoundaryFlag))
            testCase.verifyEqual(defaultSpec.variable, internalEVP.formulation)
            testCase.verifyFalse(isfield(defaultSpec, removedBoundaryFlag))
            testCase.verifyEqual(gSpec.interiorWeight(z, context), N2(z)/g, RelTol=1e-12)
            testCase.verifyEqual(fSpec.interiorWeight(z, context), ones(size(z)), AbsTol=0)
        end

        function solverReturnsInternalModesBasisWithFAndG(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);

            basisSet = solver.solveEVP(evp, nModes=3);
            z = linspace(zDomain(1), zDomain(2), 16).';
            coefficients = ones(3,1);

            testCase.verifyClass(basisSet, "IMInternalModesBasis")
            testCase.verifyEqual(basisSet.h, evp.hFromEigenvalue(basisSet.eigenvalues), RelTol=1e-12)
            testCase.verifySize(basisSet.G(z), [16 3])
            testCase.verifySize(basisSet.F(z), [16 3])
            testCase.verifyEqual(basisSet.gramMatrix(), basisSet.gramMatrix("G"), RelTol=1e-12)
            testCase.verifyEqual(size(basisSet.gramMatrix("G")), [3 3])
            testCase.verifyEqual(size(basisSet.gramMatrix("F")), [3 3])
            testCase.verifyEqual(basisSet.partialGramMatrix(zDomain(1), zDomain(2)), ...
                basisSet.partialGramMatrix("G", zDomain(1), zDomain(2)), RelTol=1e-12)
            testCase.verifyEqual(size(basisSet.partialWindowModes("F", zDomain(1), zDomain(2)).gramMatrix), [3 3])
            testCase.verifySize(basisSet.spectrum(coefficients), [3 1])
            testCase.verifySize(basisSet.spectrum(coefficients, variable="F"), [3 1])
            testCase.verifySize(basisSet.crossSpectrum(coefficients, coefficients, variable="G"), [3 1])
        end

        function internalModeCustomDepthMappingIsHonored(testCase)
            [N2, zDomain, nEVP, ~, g] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            hFromEigenvalue = @(lambda) 7 ./ lambda;
            evp = IMInternalModes(name="customDepthMap", formulation="G", N2=N2, zDomain=zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), r=@(z,ctx) ctx.N2(z)/ctx.g, ...
                g=g, hFromEigenvalue=hFromEigenvalue);

            basisSet = solver.solveEVP(evp, nModes=2);

            testCase.verifyClass(basisSet, "IMInternalModesBasis")
            testCase.verifyEqual(basisSet.h, hFromEigenvalue(basisSet.eigenvalues), RelTol=1e-12)
        end

        function evpContextOwnsDomainAndInternalModeMedium(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, f0=f0, g=g);
            configuredSolver = solver.configuredForEVP(evp);

            context = evp.contextForSolver(configuredSolver);

            testCase.verifyEqual(context.zDomain, zDomain)
            testCase.verifyEqual(context.N2(zDomain(:)), N2(zDomain(:)), RelTol=1e-12)
            testCase.verifyTrue(isfield(context, "dzLogN2"))
            testCase.verifyEqual(context.f0, f0, AbsTol=0)
            testCase.verifyEqual(context.g, g, AbsTol=0)
            testCase.verifyEqual(context.formulation, "G")
            testCase.verifyEqual(evp.parameters.f0, f0, AbsTol=0)
            testCase.verifyEqual(evp.parameters.g, g, AbsTol=0)
            testCase.verifyEqual(evp.parameters.formulation, "G")
        end

        function basisSetRetainsEVPDomainAndMedium(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);

            basisSet = solver.solveEVP(evp, nModes=2);

            testCase.verifyEqual(basisSet.zDomain, zDomain)
            testCase.verifyEqual(basisSet.N2(zDomain(:)), N2(zDomain(:)), RelTol=1e-12)
            testCase.verifyEqual(basisSet.metadata, struct())
            testCase.verifyTrue(isfield(basisSet.evp.parameters, "formulation"))
        end

        function spectralSolverSupportsCoordinateKinds(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
            coordinateKinds = ["z", "wkb", "density"];
            for coordinateKind = coordinateKinds
                solver = IMSolverSpectral(nEVP=nEVP, coordinateKind=coordinateKind);
                configuredSolver = solver.configuredForEVP(evp);
                testCase.verifyEqual(solver.coordinateKind, coordinateKind)
                testCase.verifySize(configuredSolver.physicalDerivativeMatrix(1), [nEVP nEVP])
            end
        end

        function spectralSolverRejectsUnknownCoordinateKind(testCase)
            didThrow = false;
            try
                IMSolverSpectral(coordinateKind="bad");
            catch
                didThrow = true;
            end

            testCase.verifyTrue(didThrow)
        end

        function genericEVPRejectsStratificationCoordinates(testCase)
            [~, zDomain, nEVP] = testCase.profile();
            evp = IMEigenvalueProblem(zDomain=zDomain, p=1, q=0, r=1);
            zSolver = IMSolverSpectral(nEVP=nEVP, coordinateKind="z");
            wkbSolver = IMSolverSpectral(nEVP=nEVP, coordinateKind="wkb");
            densitySolver = IMSolverSpectral(nEVP=nEVP, coordinateKind="density");

            basisSet = zSolver.solveEVP(evp, nModes=2);

            testCase.verifyFalse(isprop(evp, "hFromEigenvalue"))
            testCase.verifyClass(basisSet, "IMBasisSet")
            testCase.verifyFalse(isprop(basisSet, "h"))
            testCase.verifySize(basisSet.u(linspace(zDomain(1),zDomain(2),8).'), [8 2])
            testCase.verifySize(basisSet.uz(linspace(zDomain(1),zDomain(2),8).'), [8 2])
            testCase.verifyEqual(size(basisSet.gramMatrix()), [2 2])
            testCase.verifyEqual(size(basisSet.partialGramMatrix(zDomain(1), zDomain(2))), [2 2])
            testCase.verifyEqual(size(basisSet.partialWindowModes(zDomain(1), zDomain(2)).gramMatrix), [2 2])
            testCase.verifySize(basisSet.spectrum(ones(2,1)), [2 1])
            testCase.verifySize(basisSet.crossSpectrum(ones(2,1), ones(2,1)), [2 1])
            testCase.verifyError(@() wkbSolver.solveEVP(evp, nModes=2), "IMEigenvalueProblem:UnsupportedCoordinateKind")
            testCase.verifyError(@() densitySolver.solveEVP(evp, nModes=2), "IMEigenvalueProblem:UnsupportedCoordinateKind")
        end

        function finiteDifferenceSolverSolvesCanonicalProblem(testCase)
            z = linspace(-1000, 0, 32).';
            solver = IMSolverFiniteDifference(z=z);
            evp = IMEigenvalueProblem(zDomain=[-1000 0], p=1, q=0, r=1);

            basisSet = solver.solveEVP(evp, nModes=2);

            testCase.verifyClass(basisSet, "IMBasisSet")
            testCase.verifyFalse(isprop(basisSet, "h"))
            testCase.verifySize(basisSet.u(linspace(-1000,0,8).'), [8 2])
        end

        function finiteDifferenceSolverRejectsDomainMismatch(testCase)
            solver = IMSolverFiniteDifference(z=linspace(-1000,0,32).');
            evp = IMEigenvalueProblem(zDomain=[-900 0], p=1, q=0, r=1);

            testCase.verifyError(@() solver.solveEVP(evp, nModes=2), ...
                "IMSolverFiniteDifference:DomainMismatch")
        end

        function noValidEigenvalueDiagnosticStillReportsMatrixStats(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMEigenvalueProblem(name="degenerate", zDomain=zDomain, p=0, q=0, r=0);

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
