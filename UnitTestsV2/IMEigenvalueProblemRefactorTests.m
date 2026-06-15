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

        function internalModeValidationUsesBuiltInArgumentValidation(testCase)
            [N2, zDomain] = testCase.profile();
            badFormulationThrows = false;
            badVariableThrows = false;
            badWavenumberThrows = false;

            evp = IMInternalModes(name="charFormulation", formulation='F', N2=N2, zDomain=zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), r=@(z,~) ones(size(z)));
            try
                IMInternalModes(formulation="u", N2=N2, zDomain=zDomain);
            catch
                badFormulationThrows = true;
            end
            try
                evp.innerProduct("u");
            catch
                badVariableThrows = true;
            end
            try
                IMInternalModes.waveModesAtWavenumber(N2=N2, zDomain=zDomain, k=Inf);
            catch
                badWavenumberThrows = true;
            end

            testCase.verifyEqual(evp.formulation, "F")
            testCase.verifyTrue(badFormulationThrows)
            testCase.verifyTrue(badVariableThrows)
            testCase.verifyTrue(badWavenumberThrows)
        end

        function summarizePrintsReadableCanonicalProblem(testCase)
            [~, zDomain] = testCase.profile();
            evp = IMEigenvalueProblem(name="dirichletSummary", zDomain=zDomain, p=1, q=0, r=1, ...
                parameters=struct("alpha",2));

            output = string(evalc('evp.summarize();'));

            testCase.verifyTrue(contains(output, "dirichletSummary"))
            testCase.verifyTrue(contains(output, "-(d/dz)(p(z) du/dz) + q(z) u = lambda r(z) u"))
            testCase.verifyTrue(contains(output, "z in [-1000, 0]"))
            testCase.verifyTrue(contains(output, "surface: u(surface) = 0"))
            testCase.verifyTrue(contains(output, "bottom: u(bottom) = 0"))
            testCase.verifyTrue(contains(output, "Endpoint norm weights"))
            testCase.verifyTrue(contains(output, "none"))
            testCase.verifyTrue(contains(output, "names: alpha"))
            testCase.verifyFalse(contains(output, "function_handle"))
            testCase.verifyFalse(contains(output, "Internal-mode context"))
        end

        function summarizePrintsReadableRobinNeumannAndEndpointWeights(testCase)
            [~, zDomain] = testCase.profile();
            surface = IMBoundaryCondition(a=2, b=3);
            bottom = IMBoundaryCondition.neumann();
            activeSurface = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            robinEVP = IMEigenvalueProblem(name="robinSummary", zDomain=zDomain, p=1, q=0, r=1, ...
                surfaceBoundary=surface, bottomBoundary=bottom);
            activeEVP = IMEigenvalueProblem(name="activeSummary", zDomain=zDomain, p=1, q=0, r=1, ...
                surfaceBoundary=activeSurface);

            robinOutput = string(evalc('robinEVP.summarize();'));
            activeOutput = string(evalc('activeEVP.summarize();'));

            testCase.verifyTrue(contains(robinOutput, "surface: 2*u(surface) - 3*p(surface)*du/dz(surface) = 0"))
            testCase.verifyTrue(contains(robinOutput, "bottom: p(bottom)*du/dz(bottom) = 0"))
            testCase.verifyTrue(contains(activeOutput, "surface: -[-p(surface)*du/dz(surface)] = lambda*[u(surface)]"))
            testCase.verifyTrue(contains(activeOutput, "surface: +1 * (u(surface))^2"))
        end

        function summarizeWithSolverPrintsGridAssessment(testCase)
            [~, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            surface = IMBoundaryCondition(a=1, b=1);
            bottom = IMBoundaryCondition(a=-1, b=1);
            evp = IMEigenvalueProblem(name="gridSummary", zDomain=zDomain, p=1, q=0, r=1, ...
                surfaceBoundary=surface, bottomBoundary=bottom);

            output = string(evalc('evp.summarize(solver);'));

            testCase.verifyTrue(contains(output, "Solver assessment"))
            testCase.verifyTrue(contains(output, "solver: IMSolverSpectral"))
            testCase.verifyTrue(contains(output, "coordinate: z"))
            testCase.verifyTrue(contains(output, "grid size: 24"))
            testCase.verifyTrue(contains(output, "Coefficient ranges on solver grid"))
            testCase.verifyTrue(contains(output, "p: [1, 1]"))
            testCase.verifyTrue(contains(output, "q: [0, 0]"))
            testCase.verifyTrue(contains(output, "r: [1, 1]"))
            testCase.verifyTrue(contains(output, "negative modes: possible range [0, 2]"))
            testCase.verifyTrue(contains(output, "zero mode: absent"))
            testCase.verifyTrue(contains(output, "reason:"))
        end

        function summarizeInternalModesAppendsPhysicalContext(testCase)
            [N2, zDomain, ~, f0, g] = testCase.profile();
            evp = IMInternalModes.waveModesAtWavenumber(N2=N2, zDomain=zDomain, k=1e-4, f0=f0, g=g);

            output = string(evalc('evp.summarize();'));

            testCase.verifyTrue(contains(output, "Internal-mode context"))
            testCase.verifyTrue(contains(output, "formulation: G"))
            testCase.verifyTrue(contains(output, "f0: 0.0001 s^-1"))
            testCase.verifyTrue(contains(output, "g: 9.81 m s^-2"))
            testCase.verifyTrue(contains(output, "equivalent depth: h = hFromEigenvalue(lambda)"))
            testCase.verifyTrue(contains(output, "factory parameters: k"))
            testCase.verifyFalse(contains(output, "function_handle"))
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
            testCase.verifyFalse(isprop(basisSet, "index"))
            testCase.verifyTrue(isprop(basisSet, "modeSelectionDiagnostics"))
            testCase.verifyEqual(basisSet.modeSelectionDiagnostics.zeroModeStatus, "absent")
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
            testCase.verifyFalse(isfield(context, "dz" + "LogN2"))
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
            testCase.verifyClass(basisSet, "IMInternalModesBasis")
            testCase.verifyTrue(isprop(basisSet, "N2"))
            testCase.verifyEqual(basisSet.N2(zDomain(:)), N2(zDomain(:)), RelTol=1e-12)
            testCase.verifyFalse(ismethod(basisSet, "dz" + "LogN2"))
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
            testCase.verifyFalse(isprop(basisSet, "index"))
            testCase.verifyTrue(isprop(basisSet, "modeSelectionDiagnostics"))
            testCase.verifyFalse(isprop(basisSet, "N2"))
            testCase.verifyFalse(isprop(basisSet, "N2" + "Function"))
            testCase.verifyFalse(ismethod(basisSet, "dz" + "LogN2"))
            testCase.verifyFalse(ismethod(basisSet, "eval" + "uate"))
            testCase.verifyFalse(ismethod(basisSet, "eval" + "uateAll"))
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

        function genericNormalizationRulesUseStringNames(testCase)
            [~, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            constantScale = 3;
            evp = IMEigenvalueProblem(zDomain=zDomain, p=1, q=0, r=1);
            basisSet = solver.solveEVP(evp, nModes=2);
            basisSet = basisSet.addNormalization("constantScaled", ...
                @(basisSet,iMode) constantScale*basisSet.innerProductNormFactor(iMode));
            basisSet = basisSet.addNormalization("eigenvalueScaled", ...
                @(basisSet,iMode) sqrt(abs(basisSet.eigenvalues(iMode))) * basisSet.innerProductNormFactor(iMode));
            basisSet.normalization = "constantScaled";

            z = linspace(zDomain(1), zDomain(2), 8).';
            unityFactors = basisSet.normalizationFactors("unity");
            constantFactors = basisSet.normalizationFactors("constantScaled");
            eigenvalueFactors = basisSet.normalizationFactors("eigenvalueScaled");
            names = basisSet.normalizationNames();

            testCase.verifyFalse(isprop(evp, "normalizations"))
            testCase.verifyFalse(isprop(evp, "normalization" + "Rules"))
            testCase.verifyFalse(isprop(evp, "default" + "Normalization"))
            testCase.verifyTrue(ismember("constantScaled", names))
            testCase.verifyTrue(ismember("eigenvalueScaled", names))
            testCase.verifyEqual(basisSet.normalization, "constantScaled")
            testCase.verifyEqual(constantFactors, constantScale*unityFactors, RelTol=1e-12)
            testCase.verifyEqual(eigenvalueFactors, sqrt(abs(basisSet.eigenvalues)).*unityFactors, RelTol=1e-12)
            testCase.verifyEqual(basisSet.u(z), basisSet.u(z, normalization="unity")/constantScale, RelTol=1e-12)

            overwriteScale = 4;
            basisSet = basisSet.addNormalization("constantScaled", ...
                @(basisSet,iMode) overwriteScale*basisSet.innerProductNormFactor(iMode));
            testCase.verifyEqual(basisSet.normalizationFactors("constantScaled"), ...
                overwriteScale*unityFactors, RelTol=1e-12)
        end

        function internalModeBasisAcceptsEnumsAndAddNormalization(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes(name="customInternalNormalization", formulation="G", N2=N2, zDomain=zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), r=@(z,ctx) ctx.N2(z)/ctx.g);
            basisSet = solver.solveEVP(evp, nModes=2);
            basisSet = basisSet.addNormalization("customG", ...
                @(basisSet,iMode) 2*basisSet.innerProductNormFactor(iMode, variable="G"));
            names = basisSet.normalizationNames();

            testCase.verifyEqual(basisSet.normalization, "unity")
            testCase.verifyTrue(ismember("customG", names))
            testCase.verifyFalse(ismember("geostrophic", names))
            testCase.verifyEqual(basisSet.normalizationFactors("customG"), ...
                2*basisSet.normalizationFactors(Normalization.unity), RelTol=1e-12)
            basisSet.normalization = Normalization.unity;
            testCase.verifyEqual(basisSet.normalizationFactors(), ...
                basisSet.normalizationFactors("unity"), RelTol=1e-12)
        end

        function standardInternalModeFactoryDefaultsAreInstalledOnBasisSets(testCase)
            [N2, zDomain, nEVP, f0] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);

            hydrostaticBasis = solver.solveEVP(IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain), nModes=2);
            wavenumberBasis = solver.solveEVP(IMInternalModes.waveModesAtWavenumber(N2=N2, zDomain=zDomain, k=1e-4, f0=f0), nModes=2);
            frequencyBasis = solver.solveEVP(IMInternalModes.waveModesAtFrequency(N2=N2, zDomain=zDomain, omega=1e-3, f0=f0), nModes=2);

            testCase.verifyEqual(hydrostaticBasis.normalization, "geostrophic")
            testCase.verifyTrue(ismember("geostrophic", hydrostaticBasis.normalizationNames()))
            testCase.verifyEqual(wavenumberBasis.normalization, "kConstant")
            testCase.verifyTrue(ismember("kConstant", wavenumberBasis.normalizationNames()))
            testCase.verifyFalse(ismember("geostrophic", wavenumberBasis.normalizationNames()))
            testCase.verifyEqual(frequencyBasis.normalization, "omegaConstant")
            testCase.verifyTrue(ismember("omegaConstant", frequencyBasis.normalizationNames()))
            testCase.verifyFalse(ismember("geostrophic", frequencyBasis.normalizationNames()))
            testCase.verifyError(@() wavenumberBasis.normalizationFactors(Normalization.geostrophic), ...
                "IMBasisSet:UnsupportedNormalization")
        end

        function basisSetValidationRejectsBadScalarAnalysisInputs(testCase)
            [~, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMEigenvalueProblem(zDomain=zDomain, p=1, q=0, r=1);
            basisSet = solver.solveEVP(evp, nModes=2);
            badCoordinateThrows = false;

            try
                basisSet.u([zDomain(1); NaN]);
            catch
                badCoordinateThrows = true;
            end

            testCase.verifyTrue(badCoordinateThrows)
            testCase.verifyError(@() basisSet.partialGramMatrix(zDomain(2), zDomain(1)), "IMBasisSet:InvalidInterval")
            testCase.verifyError(@() basisSet.spectrum(ones(1,1)), "IMBasisSet:InvalidCoefficientCount")
            testCase.verifyError(@() basisSet.crossSpectrum(ones(2,1), ones(1,1)), "IMBasisSet:InvalidCoefficientCount")
            testCase.verifyEqual(basisSet.normalizationFactors(), ...
                basisSet.normalizationFactors(basisSet.normalization), RelTol=1e-12)
            testCase.verifyEqual(basisSet.normalizedNativeModes(), ...
                basisSet.normalizedNativeModes(basisSet.normalization), RelTol=1e-12)
        end

        function internalModeBasisUsesExplicitEvaluationMethods(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
            basisSet = solver.solveEVP(evp, nModes=2);

            testCase.verifyFalse(ismethod(basisSet, "eval" + "uate"))
            testCase.verifyFalse(ismethod(basisSet, "eval" + "uateAll"))
            testCase.verifySize(basisSet.F(linspace(zDomain(1),zDomain(2),8).'), [8 2])
            testCase.verifySize(basisSet.G(linspace(zDomain(1),zDomain(2),8).'), [8 2])
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
