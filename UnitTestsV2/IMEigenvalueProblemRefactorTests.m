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
        function spectralWaveAtWavenumberAssemblyMatchesPhysicalStrongForm(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            k = 1e-4;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=f0, g=g);

            [A, B] = evp.assemble(solver);
            expectedA = solver.Txx - k*k*solver.T;
            expectedB = diag((f0*f0 - N2(solver.zNative))/g)*solver.T;

            interiorRows = 2:(nEVP-1);
            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(A(1,:), solver.T(1,:), AbsTol=1e-12)
            testCase.verifyEqual(B(1,:), zeros(1,nEVP), AbsTol=1e-12)
            testCase.verifyEqual(A(end,:), solver.T(end,:), AbsTol=1e-12)
            testCase.verifyEqual(B(end,:), zeros(1,nEVP), AbsTol=1e-12)
        end

        function spectralWaveAtFrequencyAssemblyMatchesPhysicalStrongForm(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            omega = 0.8*5.2e-3;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtFrequency(omega=omega, f0=f0, g=g);

            [A, B] = evp.assemble(solver);
            expectedA = solver.Txx;
            expectedB = diag((omega*omega - N2(solver.zNative))/g)*solver.T;

            interiorRows = 2:(nEVP-1);
            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(A(1,:), solver.T(1,:), AbsTol=1e-12)
            testCase.verifyEqual(B(1,:), zeros(1,nEVP), AbsTol=1e-12)
            testCase.verifyEqual(A(end,:), solver.T(end,:), AbsTol=1e-12)
            testCase.verifyEqual(B(end,:), zeros(1,nEVP), AbsTol=1e-12)
        end

        function hydrostaticGModesDeclareHydrostaticStrongForm(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g);

            [A, B] = evp.assemble(solver);
            expectedA = solver.Txx;
            expectedB = diag(-N2(solver.zNative)/g)*solver.T;

            interiorRows = 2:(nEVP-1);
            testCase.verifyEqual(evp.name, "hydrostaticGModes")
            testCase.verifyEmpty(fieldnames(evp.parameters))
            testCase.verifyEqual(evp.defaultNormalization, Normalization.geostrophic)
            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-12)
            testCase.verifyTrue(isfield(evp.normalizations, "geostrophic"))
        end

        function waveModeMetadataSeparatesInnerWeightsBoundariesAndNormalizations(testCase)
            [~, ~, ~, f0, g] = testCase.profile();
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            legacyPrimary = "primary" + "Component";
            legacyComponents = "com" + "ponents";
            legacyEigenvalueName = "eigenvalue" + "Name";
            legacyOrdering = "order" + "ing";

            testCase.verifyEqual(evp.formulation, "G")
            testCase.verifyFalse(isprop(evp, legacyPrimary))
            testCase.verifyFalse(isprop(evp, legacyComponents))
            testCase.verifyFalse(isprop(evp, legacyEigenvalueName))
            testCase.verifyFalse(isprop(evp, legacyOrdering))
            testCase.verifyEqual(evp.hFromEigenvalue([1 2 4]), [1 0.5 0.25], AbsTol=0)
            testCase.verifyTrue(isfield(evp.innerWeights, "G"))
            testCase.verifyTrue(isfield(evp.innerWeights, "F"))
            testCase.verifyTrue(isa(evp.innerWeights.G, "function_handle"))
            testCase.verifyClass(evp.surfaceBoundary, "IMBoundary")
            testCase.verifyClass(evp.bottomBoundary, "IMBoundary")
            testCase.verifyEqual(evp.surfaceBoundary.family, "rigid")
            testCase.verifyEqual(evp.bottomBoundary.family, "rigid")
            testCase.verifyEmpty(evp.surfaceWeights)
            testCase.verifyEmpty(evp.bottomWeights)
            testCase.verifyTrue(isfield(evp.normalizations, "unity"))
            testCase.verifyTrue(isfield(evp.normalizations, "kConstant"))
            testCase.verifyTrue(isfield(evp.normalizations, "surfacePressure"))
            testCase.verifyTrue(isa(evp.normalizations.unity, "function_handle"))

            spec = evp.innerProduct("G");
            testCase.verifyEqual(spec.variable, "G")
            testCase.verifyTrue(isa(spec.interiorWeight, "function_handle"))
            testCase.verifyEmpty(spec.surfaceWeights)
            testCase.verifyEmpty(spec.bottomWeights)
            testCase.verifyTrue(spec.hasKnownBoundaryWeights)
        end

        function evpFactoriesDeclareSurfacePressureNormalization(testCase)
            [~, ~, ~, f0, g] = testCase.profile();
            evps = {
                IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g)
                IMEigenvalueProblem.waveModesAtFrequency(omega=1e-3, f0=f0, g=g)
                IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g)
                IMEigenvalueProblem.hydrostaticFModes(g=g)
            };

            for iEVP = 1:length(evps)
                evp = evps{iEVP};
                testCase.verifyTrue(isfield(evp.normalizations, "surfacePressure"))
                testCase.verifyTrue(isa(evp.normalizations.surfacePressure, "function_handle"))
            end
        end

        function defaultOperatorSyntaxCreatesStrongFormOperator(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            operator = IMOperator().plus(coefficient=@(z,ctx) ctx.N2(z), derivativeOrder=0);

            matrix = operator.matrix(solver);
            expected = diag(N2(solver.zNative))*solver.T;

            testCase.verifyEqual(operator.form, "strong")
            testCase.verifyEqual(matrix, expected, AbsTol=1e-12)
        end

        function boundaryLawsResolveImpliedVariables(testCase)
            rigidG = IMBoundary.rigid().at("surface", formulation="G");
            rigidF = IMBoundary.rigid().at("bottom", formulation="F");
            noSlipG = IMBoundary.noSlip().at("surface", formulation="G");
            noSlipF = IMBoundary.noSlip().at("bottom", formulation="F");
            dirichletF = IMBoundary.dirichlet().at("surface", formulation="F");
            neumannG = IMBoundary.neumann().at("bottom", formulation="G");
            customBoundary = IMBoundary.custom(left=IMOperator().plus(derivativeOrder=1), variable="F");
            customF = customBoundary.at("surface", formulation="G");

            testCase.verifyEqual(rigidG.variable, "G")
            testCase.verifyEqual(rigidG.family, "rigid")
            testCase.verifyEqual(rigidG.location, "surface")
            testCase.verifyEqual(rigidF.variable, "F")
            testCase.verifyEqual(rigidF.leftOperator.terms(1).derivativeOrder, 1)
            testCase.verifyEqual(noSlipG.variable, "G")
            testCase.verifyEqual(noSlipG.leftOperator.terms(1).derivativeOrder, 1)
            testCase.verifyEqual(noSlipF.variable, "F")
            testCase.verifyEqual(noSlipF.leftOperator.terms(1).derivativeOrder, 0)
            testCase.verifyEqual(dirichletF.variable, "F")
            testCase.verifyEqual(neumannG.variable, "G")
            testCase.verifyEqual(customF.family, "custom")
            testCase.verifyEqual(customF.variable, "F")
            testCase.verifyEqual(customF.leftOperator.terms(1).derivativeOrder, 1)
        end

        function evpFactoriesAcceptUpperAndLowerBoundaryLaws(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g, ...
                surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.rigid());
            evpFrequency = IMEigenvalueProblem.waveModesAtFrequency(omega=0.8*5.2e-3, f0=f0, g=g, ...
                surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.noSlip());
            evpG = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g, ...
                surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.rigid());
            evpF = IMEigenvalueProblem.hydrostaticFModes(g=g, ...
                surfaceBoundary=IMBoundary.noSlip(), bottomBoundary=IMBoundary.rigid());

            [A, B] = evp.assemble(solver);
            Dz = solver.physicalDerivativeMatrix(1);

            testCase.verifyEqual(evpFrequency.bottomBoundary.family, "noSlip")
            testCase.verifyEqual(evpFrequency.surfaceBoundary.family, "free")
            testCase.verifyEqual(evpG.bottomBoundary.family, "rigid")
            testCase.verifyEqual(evpG.surfaceBoundary.family, "free")
            testCase.verifyEqual(evpF.bottomBoundary.family, "rigid")
            testCase.verifyEqual(evpF.surfaceBoundary.family, "noSlip")
            testCase.verifyEqual(evp.bottomBoundary.family, "rigid")
            testCase.verifyEqual(evp.surfaceBoundary.family, "free")
            testCase.verifyEqual(evp.surfaceWeights(1).location, "surface")
            testCase.verifyEqual(evp.surfaceWeights(1).innerProduct, "G")
            testCase.verifyEqual(A(1,:), Dz(1,:), AbsTol=1e-12)
            testCase.verifyEqual(B(1,:), solver.T(1,:), AbsTol=1e-12)
        end

        function freeBoundaryDeclaresRowsAndBoundaryWeights(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            surfaceCondition = IMBoundary.free().at("surface", formulation="G");
            bottomCondition = IMBoundary.free().at("bottom", formulation="G");

            A = zeros(nEVP,nEVP);
            B = zeros(nEVP,nEVP);
            [A, B] = solver.applyEndpointLaw(A, B, surfaceCondition);
            [A, B] = solver.applyEndpointLaw(A, B, bottomCondition);
            surfaceWeight = surfaceCondition.boundaryWeights(1);
            bottomWeight = bottomCondition.boundaryWeights(1);
            Dz = solver.physicalDerivativeMatrix(1);

            testCase.verifyEqual(A(1,:), Dz(1,:), AbsTol=1e-12)
            testCase.verifyEqual(B(1,:), solver.T(1,:), AbsTol=1e-12)
            testCase.verifyEqual(A(end,:), Dz(end,:), AbsTol=1e-12)
            testCase.verifyEqual(B(end,:), solver.T(end,:), AbsTol=1e-12)
            testCase.verifyEqual(surfaceWeight.innerProduct, "G")
            testCase.verifyEqual(surfaceWeight.location, "surface")
            testCase.verifyEqual(bottomWeight.location, "bottom")
            testCase.verifyEqual(surfaceWeight.coefficient, 1, AbsTol=0)
            testCase.verifyEqual(bottomWeight.coefficient, -1, AbsTol=0)
            testCase.verifyEqual(surfaceWeight.leftVariable, "G")
            testCase.verifyEqual(surfaceWeight.rightVariable, "G")
        end

        function customBoundaryWeightsPlaceAndOrientOnlyLocationFreeWeights(testCase)
            locationFreeWeight = IMBoundaryWeight(innerProduct="G", coefficient=2, leftVariable="G", rightVariable="G");
            explicitBottomWeight = IMBoundaryWeight(innerProduct="G", location="bottom", ...
                coefficient=3, leftVariable="G", rightVariable="G");
            functionWeight = IMBoundaryWeight(innerProduct="G", coefficient=@(ctx) ctx.g, ...
                leftVariable="G", rightVariable="G");
            boundary = IMBoundary.custom(left=IMOperator().plus(derivativeOrder=0), ...
                boundaryWeights=[locationFreeWeight; explicitBottomWeight; functionWeight]);
            placedBoundary = boundary.at("bottom", formulation="G");

            testCase.verifyEqual(placedBoundary.boundaryWeights(1).location, "bottom")
            testCase.verifyEqual(placedBoundary.boundaryWeights(1).coefficient, -2, AbsTol=0)
            testCase.verifyEqual(placedBoundary.boundaryWeights(2).location, "bottom")
            testCase.verifyEqual(placedBoundary.boundaryWeights(2).coefficient, 3, AbsTol=0)
            testCase.verifyEqual(placedBoundary.boundaryWeights(3).coefficient(struct("g", 9.81)), -9.81, AbsTol=0)
        end

        function boundaryWeightValidatesVariablesLocationAndCoefficient(testCase)
            weight = IMBoundaryWeight(innerProduct="G", location="surface", coefficient=@(ctx) ctx.g, ...
                leftVariable="G", leftDerivativeOrder=0, rightVariable="F", rightDerivativeOrder=1);

            testCase.verifyEqual(weight.innerProduct, "G")
            testCase.verifyEqual(weight.leftVariable, "G")
            testCase.verifyEqual(weight.rightVariable, "F")
            testCase.verifyEqual(weight.rightDerivativeOrder, 1)
            testCase.verifyEqual(weight.coefficient(struct("g", 9.81)), 9.81, AbsTol=0)
            testCase.verifyError(@() IMBoundaryWeight(innerProduct="U"), "IMBoundaryWeight:InvalidVariable")
            testCase.verifyError(@() IMBoundaryWeight(location="middle"), "IMBoundaryWeight:InvalidLocation")
            testCase.verifyError(@() weight.at("bottom"), "IMBoundaryWeight:LocationMismatch")
        end

        function namedBoundaryLawsResolveForFModeRows(testCase)
            [N2, zDomain, nEVP, ~, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            freeF = IMBoundary.free().at("surface", formulation="F");
            noSlipF = IMBoundary.noSlip().at("bottom", formulation="F");
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);
            context = evp.contextForSolver(solver);

            A = zeros(nEVP,nEVP);
            B = zeros(nEVP,nEVP);
            [A, B] = solver.applyEndpointLaw(A, B, freeF, context=context);
            [A, B] = solver.applyEndpointLaw(A, B, noSlipF, context=context);
            Dz = solver.physicalDerivativeMatrix(1);
            expectedFreeRow = solver.T(1,:) + g/N2(zDomain(2))*Dz(1,:);

            testCase.verifyEqual(A(1,:), expectedFreeRow, RelTol=1e-12)
            testCase.verifyEqual(B(1,:), zeros(1,nEVP), AbsTol=0)
            testCase.verifyEqual(A(end,:), solver.T(end,:), AbsTol=1e-12)
            testCase.verifyEqual(B(end,:), zeros(1,nEVP), AbsTol=0)
        end

        function boundaryResolutionAvoidsOldHelperChain(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            source = fileread(fullfile(repoRoot, "@IMBoundary", "IMBoundary.m"));
            oldHelperNames = ["resolveAt" "resolveRigid" "resolveNoSlip" "resolveFree" "withResolvedLocation"];

            for iName = 1:length(oldHelperNames)
                testCase.verifyFalse(contains(source, oldHelperNames(iName)))
            end
        end

        function wkbCoordinateAppliesJacobianToSecondDerivative(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverWKBSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            z = solver.zNative;
            N = sqrt(N2(z));
            Nz = gradient(sqrt(N2(solver.zReference)), solver.zReference);
            Nz = interp1(solver.zReference, Nz, z, "pchip");

            D2 = solver.physicalDerivativeMatrix(2);
            expectedD2 = diag(N.*N)*solver.Txx + diag(Nz)*solver.Tx;

            testCase.verifyLessThan(norm(D2 - expectedD2, "fro")/norm(expectedD2, "fro"), 1e-10)
        end

        function densityCoordinateAppliesJacobianToSecondDerivative(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverDensitySpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            z = solver.zNative;
            N2z = gradient(N2(solver.zReference), solver.zReference);
            N2z = interp1(solver.zReference, N2z, z, "pchip");

            D2 = solver.physicalDerivativeMatrix(2);
            expectedD2 = diag(N2(z).*N2(z))*solver.Txx + diag(N2z)*solver.Tx;

            testCase.verifyLessThan(norm(D2 - expectedD2, "fro")/norm(expectedD2, "fro"), 1e-10)
        end

        function spectralUnityNormalizationUsesChebyshevInnerProduct(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            expectedFactors = testCase.expectedSpectralFactors(solver, basisSet, evp, ones(size(solver.zNative)));

            testCase.verifyEqual(basisSet.normalizationFactors(Normalization.unity), expectedFactors, RelTol=1e-11)
        end

        function wkbUnityNormalizationUsesNativeJacobian(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverWKBSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            expectedFactors = testCase.expectedSpectralFactors(solver, basisSet, evp, sqrt(N2(solver.zNative)));

            testCase.verifyEqual(basisSet.normalizationFactors(Normalization.unity), expectedFactors, RelTol=1e-11)
        end

        function densityUnityNormalizationUsesNativeJacobian(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverDensitySpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            expectedFactors = testCase.expectedSpectralFactors(solver, basisSet, evp, N2(solver.zNative));

            testCase.verifyEqual(basisSet.normalizationFactors(Normalization.unity), expectedFactors, RelTol=1e-11)
        end

        function solveEVPReturnsBasisWithChangeableNormalization(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);
            z = linspace(zDomain(1), zDomain(2), 24).';

            basisSet.normalization = Normalization.kConstant;
            Gk = basisSet.G(z);
            basisSet.normalization = Normalization.wMax;
            Gw = basisSet.G(z);
            allValues = basisSet.evaluateAll(z);

            testCase.verifySize(Gk, [length(z) 4])
            testCase.verifySize(Gw, [length(z) 4])
            testCase.verifyEqual(basisSet.evaluate("G", z), Gw, AbsTol=0)
            testCase.verifyEqual(allValues.G, Gw, AbsTol=0)
            testCase.verifyEqual(allValues.F, basisSet.F(z), AbsTol=0)
            testCase.verifyGreaterThan(norm(Gk - Gw, "fro"), 0)
        end

        function evpFactoriesDeclareWaveModeDefaultNormalization(testCase)
            [~, ~, ~, f0, g] = testCase.profile();
            wavenumberEVP = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            frequencyEVP = IMEigenvalueProblem.waveModesAtFrequency(omega=1e-3, f0=f0, g=g);

            testCase.verifyEqual(wavenumberEVP.defaultNormalization, Normalization.kConstant)
            testCase.verifyEqual(frequencyEVP.defaultNormalization, Normalization.omegaConstant)
        end

        function evpSelectionDefaultsToPositiveBaroclinicModes(testCase)
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4);

            selection = evp.selectModes([-10; 0; 3; 1; 2], 3, struct());

            testCase.verifyEqual(selection.sortIndex, [4 5 3])
            testCase.verifyEqual(selection.modeNumber, 1:3)
            testCase.verifyEqual(selection.index.positiveCount, 3)
            testCase.verifyEqual(selection.index.expectedNegativeCount, 0)
            testCase.verifyEqual(selection.index.expectedZeroCount, 0)
        end

        function evpFactoryParametersUsePhysicalInputsOnly(testCase)
            [~, ~, ~, f0, g] = testCase.profile();
            evpK = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            evpOmega = IMEigenvalueProblem.waveModesAtFrequency(omega=1e-3, f0=f0, g=g);
            evpG = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g);
            evpF = IMEigenvalueProblem.hydrostaticFModes(g=g);

            testCase.verifyEqual(evpK.name, "waveModesAtWavenumber")
            testCase.verifyEqual(evpOmega.name, "waveModesAtFrequency")
            testCase.verifyEqual(evpG.name, "hydrostaticGModes")
            testCase.verifyEqual(evpF.name, "hydrostaticFModes")
            testCase.verifyEqual(string(fieldnames(evpK.parameters)).', "k")
            testCase.verifyEqual(string(fieldnames(evpOmega.parameters)).', "omega")
            testCase.verifyEmpty(fieldnames(evpG.parameters))
            testCase.verifyEmpty(fieldnames(evpF.parameters))
            testCase.verifyEqual(evpK.parameters.k, 1e-4, AbsTol=0)
            testCase.verifyEqual(evpOmega.parameters.omega, 1e-3, AbsTol=0)

            forbiddenFields = ["problem" + "Type" "surface" + "Boundary" "bottom" + "Boundary" "f0" "g"];
            for fieldName = forbiddenFields
                testCase.verifyFalse(isfield(evpK.parameters, char(fieldName)))
                testCase.verifyFalse(isfield(evpOmega.parameters, char(fieldName)))
                testCase.verifyFalse(isfield(evpG.parameters, char(fieldName)))
                testCase.verifyFalse(isfield(evpF.parameters, char(fieldName)))
            end
            testCase.verifyFalse(isfield(evpK.parameters, "omega"))
            testCase.verifyFalse(isfield(evpOmega.parameters, "k"))
            testCase.verifyFalse(isfield(evpF.parameters, "k"))
            testCase.verifyFalse(isfield(evpF.parameters, "omega"))
            testCase.verifyEqual(evpK.f0, f0, AbsTol=0)
            testCase.verifyEqual(evpK.g, g, AbsTol=0)

            lowLevelEVP = IMEigenvalueProblem(f0=f0, g=g, parameters=struct("f0", 3, "g", 4, "k", 5));
            testCase.verifyFalse(isfield(lowLevelEVP.parameters, "f0"))
            testCase.verifyFalse(isfield(lowLevelEVP.parameters, "g"))
            testCase.verifyEqual(lowLevelEVP.parameters.k, 5, AbsTol=0)
            testCase.verifyEqual(lowLevelEVP.f0, f0, AbsTol=0)
            testCase.verifyEqual(lowLevelEVP.g, g, AbsTol=0)
        end

        function solverAndEVPContextSeparateMediumAndPhysicalConstants(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            solverContext = solver.context();
            evpContext = evp.contextForSolver(solver);

            testCase.verifyTrue(isfield(solverContext, 'N2'))
            testCase.verifyTrue(isfield(solverContext, 'dzLogN2'))
            testCase.verifyFalse(isfield(solverContext, 'g'))
            testCase.verifyFalse(isfield(solverContext, 'f0'))
            testCase.verifyTrue(isfield(solverContext, 'zDomain'))
            testCase.verifyTrue(isfield(solverContext, 'coordinateKind'))
            testCase.verifyEqual(evpContext.f0, f0, AbsTol=0)
            testCase.verifyEqual(evpContext.g, g, AbsTol=0)
        end

        function solveEVPUsesEVPDefaultNormalizationWhenOmitted(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            evp.defaultNormalization = Normalization.wMax;

            basisSet = solver.solveEVP(evp, nModes=4);

            testCase.verifyEqual(basisSet.normalization, Normalization.wMax)
        end

        function solveEVPExplainsDegenerateAssemblyWithNoValidEigenvalues(testCase)
            [N2, zDomain, nEVP, ~, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            p = @(z,ctx) ctx.f0^2 ./ ctx.N2(z);
            pz = @(z,ctx) -p(z,ctx).*ctx.dzLogN2(z);
            left = IMOperator() ...
                .plus(coefficient=p, derivativeOrder=2) ...
                .plus(coefficient=pz, derivativeOrder=1);
            right = IMOperator().plus(coefficient=@(z,ctx) -ctx.f0^2/ctx.g, derivativeOrder=0);
            evp = IMEigenvalueProblem(name="missingF0Diagnostic", formulation="F", g=g, ...
                leftOperator=left, rightOperator=right, ...
                surfaceBoundary=IMBoundary.rigid(), bottomBoundary=IMBoundary.rigid());

            try
                solver.solveEVP(evp, nModes=3);
                testCase.verifyFail("Expected IMSolver:NoValidEigenvalues.")
            catch exception
                testCase.verifyEqual(exception.identifier, 'IMSolver:NoValidEigenvalues')
                testCase.verifyTrue(contains(exception.message, "missingF0Diagnostic"))
                testCase.verifyTrue(contains(exception.message, 'norm(B,"fro")'))
                testCase.verifyTrue(contains(exception.message, "rank(B)"))
                testCase.verifyTrue(contains(exception.message, "f0"))
            end
        end

        function explicitBasisSetNormalizationOverridesEVPDefault(testCase)
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4);
            evp.defaultNormalization = Normalization.wMax;

            basisSet = IMBasisSet(evp=evp, normalization=Normalization.uMax);

            testCase.verifyEqual(basisSet.normalization, Normalization.uMax)
        end

        function basisSetFallsBackWhenEVPHasNoDefaultNormalization(testCase)
            evp = IMEigenvalueProblem();

            basisSet = IMBasisSet(evp=evp);

            testCase.verifyEqual(basisSet.normalization, Normalization.kConstant)
        end

        function basisSetDefaultsToPositionalModeNumberWhenOmitted(testCase)
            basisSet = IMBasisSet(nativeModes=zeros(3,4), h=1:4);

            testCase.verifyEqual(basisSet.modeNumber, 1:4)
        end

        function analyticalBasisFactoryUsesEVPDefaultNormalizationWhenOmitted(testCase)
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4);
            evp.defaultNormalization = Normalization.wMax;

            basisSet = IMBasisSet.constantStratification(evp=evp, nModes=4);

            testCase.verifyEqual(basisSet.normalization, Normalization.wMax)
        end

        function waveModesDoNotDefineGeostrophicNormalization(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=3);

            testCase.verifyError(@() basisSet.normalizationFactors(Normalization.geostrophic), ...
                "IMBasisSet:UnsupportedNormalization")
        end

        function hydrostaticEVPsDefaultToGeostrophicNormalization(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            basisG = solver.solveEVP(IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g), nModes=3);
            basisF = solver.solveEVP(IMEigenvalueProblem.hydrostaticFModes(g=g), nModes=3);

            testCase.verifyEqual(basisG.normalization, Normalization.geostrophic)
            testCase.verifyEqual(basisF.normalization, Normalization.geostrophic)
            testCase.verifyEqual(basisG.modeNumber, 1:3)
            testCase.verifyEqual(basisF.modeNumber, 0:2)
        end

        function solveEVPOrientsModesWithPositiveSurfaceF(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            FSurface = basisSet.evaluate("F", zDomain(2));

            testCase.verifyTrue(all(FSurface > 0))
        end

        function signConventionFallsBackToPositiveSurfaceGWhenSurfaceFIsZero(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g, ...
                surfaceBoundary=IMBoundary.noSlip(), bottomBoundary=IMBoundary.rigid());
            basisSet = solver.solveEVP(evp, nModes=4);

            FSurface = basisSet.evaluate("F", zDomain(2));
            GSurface = basisSet.evaluate("G", zDomain(2));

            testCase.verifyLessThan(max(abs(FSurface)), 1e-10)
            testCase.verifyTrue(all(GSurface > 0))
        end

        function finiteDifferenceSolverUsesSuppliedGrid(testCase)
            [N2, zDomain, ~, f0, g] = testCase.profile();
            z = linspace(zDomain(1), zDomain(2), 18).';
            solver = IMSolverFiniteDifference(z=z, N2=N2);

            testCase.verifyEqual(solver.zNative, sort(z, "descend"))
            testCase.verifyEqual(solver.physicalDerivativeMatrix(0), eye(length(z)), AbsTol=0)
            testCase.verifyEqual(solver.boundaryIndex("surface"), 1)
            testCase.verifyEqual(solver.boundaryIndex("bottom"), length(z))
            testCase.verifyTrue(isa(solver, "IMSolver"))
            basisSet = solver.solveEVP(IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g), nModes=3);
            testCase.verifyClass(basisSet, "IMBasisSet")
            testCase.verifyNumElements(basisSet.h, 3)
        end

        function partialDepthPEBoundariesCountBoundaryModes(testCase)
            positivePolicy = IMEigenvalueProblem.partialDepthPEIndexPolicy(boundarySign="positive");
            positiveIndex = positivePolicy.classify([1; 2; 3], struct());
            positiveSelection = positivePolicy.selectModes([1; 2; 3], 3, struct());

            negativePolicy = IMEigenvalueProblem.partialDepthPEIndexPolicy(boundarySign="negative");
            negativeIndex = negativePolicy.classify([-2; -1; 3], struct());
            negativeSelection = negativePolicy.selectModes([-2; -1; 3], 3, struct());

            testCase.verifyEqual(positiveIndex.expectedNegativeCount, 0)
            testCase.verifyEqual(positiveIndex.negativeCount, 0)
            testCase.verifyEqual(positiveSelection.modeNumber, [-1 -2 1])
            testCase.verifyEqual(negativeIndex.expectedNegativeCount, 2)
            testCase.verifyEqual(negativeIndex.negativeCount, 2)
            testCase.verifyEqual(negativeSelection.sortIndex, [2 1 3])
            testCase.verifyEqual(negativeSelection.modeNumber, [-1 -2 1])
        end

        function explicitBoundaryMetadataUsesEndpointModeNumbers(testCase)
            boundaryModes = [
                struct("modeNumber", -2, "indexSign", -1)
                struct("modeNumber", -1, "indexSign", 1)
            ];
            policy = IMIndexPolicy(expectedNegativeCount=1, expectedZeroCount=0, ...
                validationMode="none", boundaryModes=boundaryModes);

            selection = policy.selectModes([-1; 1; 2], 3, struct());

            testCase.verifyEqual(selection.sortIndex, [2 1 3])
            testCase.verifyEqual(selection.modeNumber, [-1 -2 1])
        end

        function customBoundaryIndexMetadataControlsBoundaryModeSelection(testCase)
            left = IMOperator().plus(derivativeOrder=0);
            boundaryMode = IMBoundary.custom(left=left, indexSign=-1, indexRank=1, boundaryModeNumber=-1);
            unlabeledIndexMetadata = IMBoundary.custom(left=left, indexSign=-1, indexRank=1);
            ordinaryBoundary = IMBoundary.custom(left=left);
            evpWithBoundaryModes = IMEigenvalueProblem(surfaceBoundary=boundaryMode, bottomBoundary=ordinaryBoundary);
            evpWithoutBoundaryModes = IMEigenvalueProblem(surfaceBoundary=unlabeledIndexMetadata, ...
                bottomBoundary=ordinaryBoundary);

            withSelection = evpWithBoundaryModes.selectModes([-10; -1; 2; 5], 3, struct());
            withoutSelection = evpWithoutBoundaryModes.selectModes([-10; -1; 2; 5], 2, struct());

            testCase.verifyEqual(withSelection.sortIndex, [2 3 4])
            testCase.verifyEqual(withSelection.modeNumber, [-1 1 2])
            testCase.verifyEqual(withSelection.index.expectedNegativeCount, 1)
            testCase.verifyEqual(withoutSelection.sortIndex, [3 4])
            testCase.verifyEqual(withoutSelection.modeNumber, 1:2)
            testCase.verifyEqual(withoutSelection.index.expectedNegativeCount, 0)
        end

        function indexPolicyErrorsWhenObservedIndexDisagrees(testCase)
            policy = IMEigenvalueProblem.partialDepthPEIndexPolicy(boundarySign="negative");

            testCase.verifyError(@() policy.classify([-1; 2; 3], struct()), ...
                "IMIndexPolicy:IndexMismatch")
        end

        function indexPolicySelectsAndLabelsNamedModeBlocks(testCase)
            policy = IMIndexPolicy.fixed(expectedNegativeCount=2, expectedZeroCount=1, validationMode="none");

            selection = policy.selectModes([-10; -1; 0; 2; 5], 5, struct());

            testCase.verifyEqual(selection.sortIndex, [2 1 3 4 5])
            testCase.verifyEqual(selection.modeNumber, [-1 -2 0 1 2])
            testCase.verifyEqual(selection.index.negativeCount, 2)
            testCase.verifyEqual(selection.index.zeroCount, 1)
        end

        function indexPolicyExcludesUnexpectedNegativeArtifacts(testCase)
            policy = IMIndexPolicy.fixed(expectedNegativeCount=0, expectedZeroCount=1, validationMode="none");

            selection = policy.selectModes([-1e20; 0; 1; 2], 3, struct());

            testCase.verifyEqual(selection.sortIndex, [2 3 4])
            testCase.verifyEqual(selection.modeNumber, [0 1 2])
            testCase.verifyEqual(selection.index.negativeCount, 0)
            testCase.verifyEqual(selection.index.zeroCount, 1)
        end

        function hydrostaticFModeEVPDeclaresBarotropicMode(testCase)
            [~, ~, ~, ~, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);

            index = evp.classifyEigenvalues([0; 1; 2; 3], struct());

            testCase.verifyTrue(evp.hasBarotropicMode)
            testCase.verifyEqual(evp.indexValidationMode, "warning")
            testCase.verifyEqual(index.expectedZeroCount, 1)
            testCase.verifyEqual(index.zeroCount, 1)
            testCase.verifyEqual(index.positiveCount, 3)
        end

        function hydrostaticFModeEVPDoesNotDeclareWavenumber(testCase)
            [~, ~, ~, ~, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);

            testCase.verifyFalse(isfield(evp.parameters, "k"))
            testCase.verifyFalse(isfield(evp.parameters, "omega"))
        end

        function hydrostaticFModeEVPEvaluatesDiagnosticGFromFDerivative(testCase)
            [N2, zDomain, nEVP, ~, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);
            warningState = warning("off", "IMIndexPolicy:IndexMismatch");
            cleanup = onCleanup(@() warning(warningState));
            basisSet = solver.solveEVP(evp, nModes=4);
            z = linspace(zDomain(1), zDomain(2), 24).';

            Fz = solver.evaluatePhysicalDerivative(basisSet.nativeModes, z, 1);
            rawGExpected = (-g./N2(z)).*Fz;
            factors = basisSet.normalizationFactors(Normalization.unity);
            GExpected = rawGExpected ./ factors;
            legacyComponents = "com" + "ponents";

            testCase.verifyEqual(evp.formulation, "F")
            testCase.verifyFalse(isprop(evp, legacyComponents))
            testCase.verifyEqual(basisSet.G(z, normalization=Normalization.unity), GExpected, RelTol=1e-10)
            testCase.verifyEqual(basisSet.evaluate("G", z, normalization=Normalization.unity), GExpected, RelTol=1e-10)
        end

        function hydrostaticGModeGeostrophicNormalizationGivesParsevalMetrics(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            basisSet = solver.solveEVP(IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g), nModes=4);

            gramG = basisSet.gramMatrix("G");
            gramF = basisSet.gramMatrix("F");

            testCase.verifyEqual(diag(gramG).', ones(1,4), AbsTol=1e-10)
            testCase.verifyEqual(diag(gramF).', basisSet.h, RelTol=1e-6, AbsTol=1e-12)
            testCase.verifyGreaterThan(basisSet.eigenvalues, zeros(size(basisSet.eigenvalues)))
        end

        function hydrostaticFModeGeostrophicNormalizationKeepsNullMode(testCase)
            [N2, zDomain, nEVP, ~, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            warningState = warning("off", "IMIndexPolicy:IndexMismatch");
            cleanup = onCleanup(@() warning(warningState));
            basisSet = solver.solveEVP(IMEigenvalueProblem.hydrostaticFModes(g=g), nModes=5);
            z = linspace(zDomain(1), zDomain(2), 64).';
            F = basisSet.evaluate("F", z);
            G = basisSet.evaluate("G", z);
            gramF = basisSet.gramMatrix("F");
            gramG = basisSet.gramMatrix("G");
            depth = diff(zDomain);

            testCase.verifyLessThan(abs(basisSet.eigenvalues(1)), 1e-10)
            testCase.verifyEqual(basisSet.modeNumber, 0:4)
            testCase.verifyTrue(isinf(abs(basisSet.h(1))))
            testCase.verifyEqual(F(:,1), ones(size(z)), AbsTol=1e-10)
            testCase.verifyLessThan(max(abs(G(:,1))), 1e-9)
            testCase.verifyEqual(gramF(1,1), depth, RelTol=1e-10)
            testCase.verifyEqual(diag(gramF(2:end,2:end)).', basisSet.h(2:end), RelTol=1e-6, AbsTol=1e-12)
            testCase.verifyEqual(diag(gramG(2:end,2:end)).', ones(1,4), AbsTol=1e-10)
        end

        function hydrostaticFAndGBaroclinicModesMatchUnderGeostrophicNormalization(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            warningState = warning("off", "IMIndexPolicy:IndexMismatch");
            cleanup = onCleanup(@() warning(warningState));
            basisF = solver.solveEVP(IMEigenvalueProblem.hydrostaticFModes(g=g), nModes=5);
            basisG = solver.solveEVP(IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g), nModes=4);
            z = linspace(zDomain(1), zDomain(2), 128).';
            FFromF = basisF.evaluate("F", z);
            FFromG = basisG.evaluate("F", z);
            GFromF = basisF.evaluate("G", z);
            GFromG = basisG.evaluate("G", z);

            testCase.verifyEqual(basisF.h(2:end), basisG.h, RelTol=1e-6, AbsTol=1e-12)
            testCase.verifyEqual(basisF.modeNumber, 0:4)
            testCase.verifyEqual(basisG.modeNumber, 1:4)
            testCase.verifyLessThan(norm(FFromF(:,2:end) - FFromG, "fro")/norm(FFromG, "fro"), 1e-5)
            testCase.verifyLessThan(norm(GFromF(:,2:end) - GFromG, "fro")/norm(GFromG, "fro"), 1e-5)
        end

        function fullDepthPartialGramMatchesGramMatrix(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            fullGram = basisSet.gramMatrix("G");
            partialGram = basisSet.partialGramMatrix("G", zDomain(1), zDomain(2));

            testCase.verifyEqual(partialGram, fullGram, AbsTol=0)
            testCase.verifyEqual(partialGram, partialGram.', AbsTol=1e-12)
        end

        function unityNormalizationMakesSpectralGramUnitDiagonal(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            basisSet = solver.solveEVP(evp, nModes=4);

            basisSet.normalization = Normalization.unity;
            gram = basisSet.gramMatrix("G");

            testCase.verifyEqual(diag(gram).', ones(1,4), AbsTol=1e-10)
        end

        function freeBoundaryInnerProductAddsSurfaceWeightOutsideIntegral(testCase)
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            zDomain = [-1300 0];
            nEVP = 48;
            g = 9.81;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0, g=g, ...
                surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.rigid());
            basisSet = solver.solveEVP(evp, nModes=3);
            z = solver.innerProductGrid(zDomain);
            rawG = solver.T*basisSet.nativeModes;
            weight = IMOperator.evaluateCoefficient(evp.innerWeights.G, z, testCase.evpContext(solver, evp));
            interior = solver.integrateInnerProduct(z, weight.*rawG(:,1).*rawG(:,1), zDomain);
            expectedFactor = sqrt(abs(interior + rawG(1,1)*rawG(1,1)));

            factors = basisSet.normalizationFactors(Normalization.unity);
            basisSet.normalization = Normalization.unity;
            fullGram = basisSet.gramMatrix("G");
            zCut = zDomain(2) - 1;
            partialGram = basisSet.partialGramMatrix("G", zDomain(1), zCut);
            GPartial = basisSet.evaluate("G", z);
            expectedPartial = solver.integrateInnerProduct(z, weight.*GPartial(:,1).*GPartial(:,1), [zDomain(1) zCut]);

            testCase.verifyEqual(basisSet.modeNumber, [-1 1 2])
            testCase.verifyEqual(factors(1), expectedFactor, RelTol=1e-10)
            testCase.verifyEqual(fullGram(1,1), 1, AbsTol=1e-10)
            testCase.verifyEqual(partialGram(1,1), expectedPartial, AbsTol=1e-10)
            testCase.verifyGreaterThan(fullGram(1,1), partialGram(1,1))
        end

        function fixedFrequencyKConstantUsesBoundaryWeights(testCase)
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            zDomain = [-1300 0];
            nEVP = 48;
            f0 = 1e-4;
            g = 9.81;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtFrequency(omega=0.8*N0, f0=f0, g=g, ...
                surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.free());
            basisSet = solver.solveEVP(evp, nModes=3);
            z = solver.innerProductGrid(zDomain);
            rawG = solver.T*basisSet.nativeModes;
            weight = IMOperator.evaluateCoefficient(@(z,ctx) (ctx.N2(z) - ctx.f0*ctx.f0)/ctx.g, ...
                z, testCase.evpContext(solver, evp));
            interior = solver.integrateInnerProduct(z, weight.*rawG(:,1).*rawG(:,1), zDomain);
            expectedFactor = sqrt(abs(interior + rawG(1,1)*rawG(1,1) - rawG(end,1)*rawG(end,1)));

            factors = basisSet.normalizationFactors(Normalization.kConstant);

            testCase.verifyEqual(factors(1), expectedFactor, RelTol=1e-10)
        end

        function mixedBoundaryWeightsContributeToGramMatrix(testCase)
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            zDomain = [-1300 0];
            nEVP = 48;
            g = 9.81;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            mixedWeights = [
                IMBoundaryWeight(innerProduct="G", location="surface", coefficient=1, ...
                    leftVariable="G", rightVariable="G")
                IMBoundaryWeight(innerProduct="G", location="surface", coefficient=1, ...
                    leftVariable="G", rightVariable="F")
                IMBoundaryWeight(innerProduct="G", location="surface", coefficient=1, ...
                    leftVariable="F", rightVariable="G")
            ];
            freeSurface = IMBoundary.custom( ...
                left=IMOperator().plus(derivativeOrder=1), ...
                right=IMOperator().plus(derivativeOrder=0), ...
                boundaryWeights=mixedWeights, ...
                indexSign=1, indexRank=1, boundaryModeNumber=-1);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0, g=g, ...
                surfaceBoundary=freeSurface, bottomBoundary=IMBoundary.rigid());
            basisSet = solver.solveEVP(evp, nModes=3);
            basisSet.normalization = Normalization.unity;

            z = solver.innerProductGrid(zDomain);
            weight = IMOperator.evaluateCoefficient(evp.innerWeights.G, z, testCase.evpContext(solver, evp));
            G = basisSet.evaluate("G", z);
            interior = solver.integrateInnerProduct(z, weight.*G(:,1).*G(:,1), zDomain);
            GSurface = basisSet.evaluate("G", zDomain(2));
            FSurface = basisSet.evaluate("F", zDomain(2));
            expected = interior + GSurface(1,1)*GSurface(1,1) ...
                + GSurface(1,1)*FSurface(1,1) + FSurface(1,1)*GSurface(1,1);

            gram = basisSet.gramMatrix("G");

            testCase.verifyEqual(gram(1,1), expected, AbsTol=1e-10)
        end

        function linearBoundaryFamiliesDeclareTrustedBoundaryWeights(testCase)
            linearF = IMBoundary.linearF(c=2, d=3).at("surface");
            linearG = IMBoundary.linearG(a=2, b=3).at("surface");
            bottomLinearF = IMBoundary.linearF(c=2, d=3).at("bottom");
            bottomLinearG = IMBoundary.linearG(a=2, b=3).at("bottom");

            testCase.verifyTrue(linearF.hasKnownBoundaryWeights)
            testCase.verifyEqual(linearF.boundaryWeights(1).innerProduct, "G")
            testCase.verifyEqual(linearF.boundaryWeights(1).location, "surface")
            testCase.verifyEqual(linearF.boundaryWeights(1).coefficient, -1.5, AbsTol=0)
            testCase.verifyEqual(bottomLinearF.boundaryWeights(1).location, "bottom")
            testCase.verifyEqual(bottomLinearF.boundaryWeights(1).coefficient, 1.5, AbsTol=0)
            testCase.verifyTrue(linearG.hasKnownBoundaryWeights)
            testCase.verifyEqual(linearG.boundaryWeights(1).innerProduct, "G")
            testCase.verifyEqual(linearG.boundaryWeights(1).coefficient, 1.5, AbsTol=0)
            testCase.verifyEqual(bottomLinearG.boundaryWeights(1).location, "bottom")
            testCase.verifyEqual(bottomLinearG.boundaryWeights(1).coefficient, -1.5, AbsTol=0)
        end

        function unresolvedLinearFamiliesWarnButStillAssembleRows(testCase)
            testCase.verifyWarning(@() IMBoundary.linearF(a=1, b=2, c=3).at("surface"), ...
                "IMBoundary:UnknownBoundaryWeights")
            warningState = warning("off", "IMBoundary:UnknownBoundaryWeights");
            cleanup = onCleanup(@() warning(warningState));
            boundary = IMBoundary.linearF(a=1, b=2, c=3).at("surface");

            testCase.verifyFalse(boundary.hasKnownBoundaryWeights)
            testCase.verifyEmpty(boundary.boundaryWeights)
            testCase.verifyError(@() IMBoundary.linearF(c=1).at("surface", formulation="G"), ...
                "IMBoundary:UnsupportedPlacement")
        end

        function finiteDifferenceUnityNormalizationUsesBoundedTrapz(testCase)
            [N2, zDomain, ~, f0, g] = testCase.profile();
            solver = IMSolverFiniteDifference(z=linspace(zDomain(1), zDomain(2), 35).', N2=N2);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g, ...
                surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.rigid());
            basisSet = solver.solveEVP(evp, nModes=3);
            z = solver.innerProductGrid(zDomain);
            rawG = solver.evaluateNativeModes(basisSet.nativeModes, z);
            weight = IMOperator.evaluateCoefficient(evp.innerWeights.G, z, testCase.evpContext(solver, evp));
            normValue = trapz(z, weight.*rawG(:,1).*rawG(:,1));
            surfaceValue = solver.evaluateNativeModes(basisSet.nativeModes(:,1), zDomain(2));
            expectedFactor = sqrt(abs(normValue + surfaceValue*surfaceValue));

            factors = basisSet.normalizationFactors(Normalization.unity);

            testCase.verifyEqual(factors(1), expectedFactor, RelTol=1e-10)
        end

        function surfacePressureNormalizationRequiresNonzeroSurfaceF(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g, ...
                surfaceBoundary=IMBoundary.noSlip(), bottomBoundary=IMBoundary.rigid());
            basisSet = solver.solveEVP(evp, nModes=3);

            testCase.verifyLessThan(max(abs(basisSet.F(zDomain(2), normalization=Normalization.uMax))), 1e-8)
            testCase.verifyError(@() basisSet.normalizationFactors(Normalization.surfacePressure), ...
                "IMBasisSet:UnsupportedNormalization")
        end

        function constantStratificationFactoryReturnsAnalyticalBasisSet(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            nModes = 6;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0);

            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);

            testCase.verifyClass(basisSet, "IMBasisSetConstantStratification")
            testCase.verifyNumElements(basisSet.h, nModes)
            testCase.verifyEqual(basisSet.eigenvalues, 1./basisSet.h, RelTol=1e-12)
            testCase.verifyEqual(basisSet.modeNumber, 1:nModes)
        end

        function constantStratificationInfersParametersFromEVP(testCase)
            N0 = 5.2e-3;
            f0 = 1e-4;
            g = 10;
            k = 1e-4;
            zDomain = [-1300 0];
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=f0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=2);
            D = diff(zDomain);
            hExpected = (N0*N0 - f0*f0)/(g*(k*k + (pi/D)^2));

            testCase.verifyEqual(basisSet.evp.f0, f0, AbsTol=0)
            testCase.verifyEqual(basisSet.evp.g, g, AbsTol=0)
            testCase.verifyEqual(basisSet.h(1), hExpected, RelTol=1e-12)
        end

        function constantStratificationMatchesDirectFixedWavenumberSolution(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            z = linspace(zDomain(1), zDomain(2), 65).';
            nModes = 6;
            k = 1e-4;
            g = 9.81;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);

            direct = InternalModesConstantStratification(N0=N0, zIn=zDomain, zOut=z, latitude=0, nModes=nModes, g=g);
            direct.normalization = Normalization.kConstant;
            [FExpected, GExpected, hExpected] = direct.modesAtWavenumber(k);

            testCase.verifyEqual(basisSet.h, hExpected, RelTol=1e-12)
            testCase.verifyEqual(basisSet.evaluate("G", z), GExpected, AbsTol=1e-10)
            testCase.verifyEqual(basisSet.evaluate("F", z), FExpected, AbsTol=1e-10)
        end

        function constantStratificationMatchesDirectFixedFrequencySolution(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            z = linspace(zDomain(1), zDomain(2), 65).';
            nModes = 6;
            omega = 0.8*N0;
            g = 9.81;
            evp = IMEigenvalueProblem.waveModesAtFrequency(omega=omega, f0=0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);

            direct = InternalModesConstantStratification(N0=N0, zIn=zDomain, zOut=z, latitude=0, nModes=nModes, g=g);
            direct.normalization = Normalization.kConstant;
            [FExpected, GExpected, hExpected] = direct.modesAtFrequency(omega);

            testCase.verifyEqual(basisSet.h, hExpected, RelTol=1e-12)
            testCase.verifyEqual(basisSet.evaluate("G", z, normalization=Normalization.kConstant), GExpected, AbsTol=1e-10)
            testCase.verifyEqual(basisSet.evaluate("F", z, normalization=Normalization.kConstant), FExpected, AbsTol=1e-10)
        end

        function constantStratificationFreeBoundaryMatchesDirectSolution(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            z = linspace(zDomain(1), zDomain(2), 65).';
            nModes = 6;
            k = 1e-4;
            g = 9.81;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=k, f0=0, g=g, ...
                surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.rigid());
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);

            direct = InternalModesConstantStratification(N0=N0, zIn=zDomain, zOut=z, latitude=0, nModes=nModes, g=g);
            direct.(sprintf("%sBoundary", "upper")) = UpperBoundary.freeSurface;
            direct.normalization = Normalization.kConstant;
            [FExpected, GExpected, hExpected] = direct.modesAtWavenumber(k);

            testCase.verifyEqual(basisSet.h, hExpected, RelTol=1e-10)
            testCase.verifyEqual(basisSet.modeNumber, [-1 1:(nModes-1)])
            testCase.verifyEqual(basisSet.evaluate("G", z), GExpected, AbsTol=1e-8)
            testCase.verifyEqual(basisSet.evaluate("F", z), FExpected, AbsTol=1e-8)
        end

        function constantStratificationNormalizationCanChange(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            z = linspace(zDomain(1), zDomain(2), 65).';
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=4);

            basisSet.normalization = Normalization.kConstant;
            Gk = basisSet.evaluate("G", z);
            basisSet.normalization = Normalization.wMax;
            Gw = basisSet.evaluate("G", z);

            testCase.verifyGreaterThan(norm(Gk - Gw, "fro"), 0)
        end

        function constantStratificationGramUsesEVPInnerWeight(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            nModes = 4;
            g = 9.81;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0, g=g);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);
            z = linspace(zDomain(1), zDomain(2), 256).';
            G = basisSet.evaluate("G", z);
            weight = (N0*N0/g)*ones(size(z));
            expectedGram = zeros(nModes,nModes);
            for iMode = 1:nModes
                for jMode = iMode:nModes
                    value = trapz(z, weight.*G(:,iMode).*G(:,jMode));
                    expectedGram(iMode,jMode) = value;
                    expectedGram(jMode,iMode) = value;
                end
            end

            testCase.verifyEqual(basisSet.gramMatrix("G"), expectedGram, AbsTol=1e-12)
        end

        function unityNormalizationMakesConstantStratificationGramUnitDiagonal(testCase)
            N0 = 5.2e-3;
            zDomain = [-1300 0];
            nModes = 4;
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0);
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);

            basisSet.normalization = Normalization.unity;
            gram = basisSet.gramMatrix("G");

            testCase.verifyEqual(diag(gram).', ones(1,nModes), AbsTol=1e-12)
            testCase.verifyLessThan(norm(gram - diag(diag(gram)), "fro"), 1e-12)
        end

        function deferredAnalyticalFactoriesThrowExplicitErrors(testCase)
            testCase.verifyError(@() IMBasisSet.wkbApproximation(), ...
                "IMBasisSet:AnalyticalBasisNotImplemented")
        end
    end

    methods (Access = private)
        function [N2, zDomain, nEVP, f0, g] = profile(~)
            N0 = 5.2e-3;
            b = 1300;
            N2 = @(z) N0*N0*exp(2*z/b);
            zDomain = [-4000 0];
            nEVP = 32;
            f0 = 1e-4;
            g = 9.81;
        end

        function factors = expectedSpectralFactors(testCase, solver, basisSet, evp, q)
            z = solver.zNative;
            rawG = solver.T*basisSet.nativeModes;
            weight = IMOperator.evaluateCoefficient(evp.innerWeights.G, z, testCase.evpContext(solver, evp));
            weights = testCase.chebyshevIntegrationWeights(solver);
            [surfaceWeight, bottomWeight] = testCase.scalarEndpointWeights(evp, "G");
            factors = zeros(1,size(rawG,2));
            for iMode = 1:size(rawG,2)
                integrandCheb = InternalModesSpectral.fct(weight.*rawG(:,iMode).*rawG(:,iMode)./q(:));
                normValue = sum(weights.*integrandCheb);
                normValue = normValue + surfaceWeight*rawG(1,iMode)*rawG(1,iMode);
                normValue = normValue + bottomWeight*rawG(end,iMode)*rawG(end,iMode);
                factors(iMode) = sqrt(abs(normValue));
            end
        end

        function [surfaceWeight, bottomWeight] = scalarEndpointWeights(~, evp, variable)
            surfaceWeight = 0;
            bottomWeight = 0;
            spec = evp.innerProduct(variable);
            weights = [spec.bottomWeights; spec.surfaceWeights];
            for iWeight = 1:length(weights)
                weight = weights(iWeight);
                if string(weight.innerProduct) ~= string(variable)
                    continue;
                end
                if ~isnumeric(weight.coefficient) || ~isscalar(weight.coefficient)
                    continue;
                end
                if string(weight.leftVariable) ~= string(variable) ...
                        || string(weight.rightVariable) ~= string(variable)
                    continue;
                end
                if weight.leftDerivativeOrder ~= 0 || weight.rightDerivativeOrder ~= 0
                    continue;
                end
                if string(weight.location) == "surface"
                    surfaceWeight = surfaceWeight + weight.coefficient;
                elseif string(weight.location) == "bottom"
                    bottomWeight = bottomWeight + weight.coefficient;
                end
            end
        end

        function context = evpContext(~, solver, evp)
            context = evp.contextForSolver(solver);
        end

        function weights = chebyshevIntegrationWeights(~, solver)
            np = (0:(solver.nEVP-1)).';
            weights = -(1+(-1).^np)./(np.*np - 1);
            weights(2) = 0;
            weights = (max(solver.xNative) - min(solver.xNative))*weights/2;
        end
    end
end
