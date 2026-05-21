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
            testCase.verifyEqual(evp.parameters.problemType, "hydrostaticGModes")
            testCase.verifyEqual(evp.defaultNormalization, Normalization.geostrophic)
            testCase.verifyEqual(A(interiorRows,:), expectedA(interiorRows,:), AbsTol=1e-12)
            testCase.verifyEqual(B(interiorRows,:), expectedB(interiorRows,:), AbsTol=1e-12)
            testCase.verifyTrue(isfield(evp.normalizations, "geostrophic"))
        end

        function waveModeMetadataSeparatesInnerWeightsBoundariesAndNormalizations(testCase)
            [~, ~, ~, f0, g] = testCase.profile();
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);

            testCase.verifyTrue(isfield(evp.components.G, "innerWeight"))
            testCase.verifyTrue(isa(evp.components.G.innerWeight, "function_handle"))
            testCase.verifyFalse(isfield(evp.components.G, "surfaceWeight"))
            testCase.verifyFalse(isfield(evp.components.G, "bottomWeight"))
            testCase.verifyClass(evp.boundaryRows, "IMBoundaryRow")
            testCase.verifyEqual([evp.boundaryRows.family], ["rigid" "rigid"])
            testCase.verifyFalse(isfield(evp.components.G, "gramWeight"))
            testCase.verifyFalse(isfield(evp.components.G, "normalizationWeights"))
            testCase.verifyFalse(isfield(evp.components.G, "normalizationBoundaryWeights"))
            testCase.verifyTrue(isfield(evp.normalizations, "unity"))
            testCase.verifyTrue(isfield(evp.normalizations, "kConstant"))
            testCase.verifyTrue(isa(evp.normalizations.unity, "function_handle"))
        end

        function boundaryLawsResolveImpliedComponents(testCase)
            rigidG = IMBoundary.rigid().resolve(endpoint="upper", primaryComponent="G");
            rigidF = IMBoundary.rigid().resolve(endpoint="lower", primaryComponent="F");
            noSlipG = IMBoundary.noSlip().resolve(endpoint="upper", primaryComponent="G");
            noSlipF = IMBoundary.noSlip().resolve(endpoint="lower", primaryComponent="F");
            dirichletPrimaryF = IMBoundary.dirichlet().resolve(endpoint="upper", primaryComponent="F");
            neumannPrimaryG = IMBoundary.neumann().resolve(endpoint="lower", primaryComponent="G");
            dirichletF = IMBoundary.dirichlet(on="F").resolve(endpoint="upper", primaryComponent="G");
            neumannG = IMBoundary.neumann(on="G").resolve(endpoint="lower", primaryComponent="F");
            operatorBoundary = IMBoundary.operator(left=IMOperator.strong().plus(derivativeOrder=1), on="F");
            operatorF = operatorBoundary.resolve(endpoint="upper", primaryComponent="G");

            testCase.verifyEqual(rigidG.component, "G")
            testCase.verifyEqual(rigidG.family, "rigid")
            testCase.verifyEqual(rigidG.endpoint, "upper")
            testCase.verifyEqual(rigidF.component, "F")
            testCase.verifyEqual(rigidF.leftOperator.terms(1).derivativeOrder, 1)
            testCase.verifyEqual(noSlipG.component, "G")
            testCase.verifyEqual(noSlipG.leftOperator.terms(1).derivativeOrder, 1)
            testCase.verifyEqual(noSlipF.component, "F")
            testCase.verifyEqual(noSlipF.leftOperator.terms(1).derivativeOrder, 0)
            testCase.verifyEqual(dirichletPrimaryF.component, "F")
            testCase.verifyEqual(neumannPrimaryG.component, "G")
            testCase.verifyEqual(dirichletF.component, "F")
            testCase.verifyEqual(neumannG.component, "G")
            testCase.verifyEqual(operatorF.family, "operator")
            testCase.verifyEqual(operatorF.component, "F")
            testCase.verifyEqual(operatorF.leftOperator.terms(1).derivativeOrder, 1)
        end

        function evpFactoriesAcceptUpperAndLowerBoundaryLaws(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g, ...
                upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.rigid());
            evpFrequency = IMEigenvalueProblem.waveModesAtFrequency(omega=0.8*5.2e-3, f0=f0, g=g, ...
                upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.noSlip());
            evpG = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g, ...
                upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.rigid());
            evpF = IMEigenvalueProblem.hydrostaticFModes(g=g, ...
                upperBoundary=IMBoundary.noSlip(), lowerBoundary=IMBoundary.rigid());

            [A, B] = evp.assemble(solver);
            Dz = solver.physicalDerivativeMatrix(1);

            testCase.verifyEqual(evpFrequency.parameters.upperBoundary, "free")
            testCase.verifyEqual(evpFrequency.parameters.lowerBoundary, "noSlip")
            testCase.verifyEqual(evpG.parameters.upperBoundary, "free")
            testCase.verifyEqual(evpF.parameters.upperBoundary, "noSlip")
            testCase.verifyEqual(evp.parameters.lowerBoundary, "rigid")
            testCase.verifyEqual(evp.parameters.upperBoundary, "free")
            testCase.verifyEqual([evp.boundaryRows.endpoint], ["lower" "upper"])
            testCase.verifyEqual(A(1,:), Dz(1,:), AbsTol=1e-12)
            testCase.verifyEqual(B(1,:), solver.T(1,:), AbsTol=1e-12)
        end

        function freeBoundaryDeclaresRowsAndInnerProductTerms(testCase)
            [N2, zDomain, nEVP, f0, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            upperBoundary = IMBoundary.free().resolve(endpoint="upper", primaryComponent="G");
            lowerBoundary = IMBoundary.free().resolve(endpoint="lower", primaryComponent="G");

            A = zeros(nEVP,nEVP);
            B = zeros(nEVP,nEVP);
            [A, B] = upperBoundary.apply(A, B, solver);
            [A, B] = lowerBoundary.apply(A, B, solver);
            upperTerm = upperBoundary.innerProductTerms(1);
            lowerTerm = lowerBoundary.innerProductTerms(1);
            Dz = solver.physicalDerivativeMatrix(1);

            testCase.verifyEqual(A(1,:), Dz(1,:), AbsTol=1e-12)
            testCase.verifyEqual(B(1,:), solver.T(1,:), AbsTol=1e-12)
            testCase.verifyEqual(A(end,:), Dz(end,:), AbsTol=1e-12)
            testCase.verifyEqual(B(end,:), solver.T(end,:), AbsTol=1e-12)
            testCase.verifyEqual(upperTerm.innerProductComponent, "G")
            testCase.verifyEqual(upperTerm.location, "surface")
            testCase.verifyEqual(lowerTerm.location, "bottom")
            testCase.verifyEqual(upperTerm.coefficient, 1, AbsTol=0)
            testCase.verifyEqual(upperTerm.leftTrace.component, "G")
            testCase.verifyEqual(upperTerm.rightTrace.component, "G")
        end

        function namedBoundaryLawsResolveForFModeRows(testCase)
            [N2, zDomain, nEVP, ~, g] = testCase.profile();
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            freeF = IMBoundary.free().resolve(endpoint="upper", primaryComponent="F");
            noSlipF = IMBoundary.noSlip().resolve(endpoint="lower", primaryComponent="F");
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);
            context = evp.contextForSolver(solver);

            A = zeros(nEVP,nEVP);
            B = zeros(nEVP,nEVP);
            [A, B] = freeF.apply(A, B, solver, context=context);
            [A, B] = noSlipF.apply(A, B, solver, context=context);
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
            Gk = basisSet.evaluate("G", z);
            basisSet.normalization = Normalization.wMax;
            Gw = basisSet.evaluate("G", z);

            testCase.verifySize(Gk, [length(z) 4])
            testCase.verifySize(Gw, [length(z) 4])
            testCase.verifyGreaterThan(norm(Gk - Gw, "fro"), 0)
        end

        function evpFactoriesDeclareWaveModeDefaultNormalization(testCase)
            [~, ~, ~, f0, g] = testCase.profile();
            wavenumberEVP = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            frequencyEVP = IMEigenvalueProblem.waveModesAtFrequency(omega=1e-3, f0=f0, g=g);

            testCase.verifyEqual(wavenumberEVP.defaultNormalization, Normalization.kConstant)
            testCase.verifyEqual(frequencyEVP.defaultNormalization, Normalization.omegaConstant)
        end

        function evpFactoryParametersUseDocumentedMetadataSchema(testCase)
            [~, ~, ~, f0, g] = testCase.profile();
            evpK = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g);
            evpOmega = IMEigenvalueProblem.waveModesAtFrequency(omega=1e-3, f0=f0, g=g);
            evpG = IMEigenvalueProblem.hydrostaticGModes(f0=f0, g=g);
            evpF = IMEigenvalueProblem.hydrostaticFModes(g=g);

            requiredFields = ["problemType" "upperBoundary" "lowerBoundary"];
            for fieldName = requiredFields
                testCase.verifyTrue(isfield(evpK.parameters, char(fieldName)))
                testCase.verifyTrue(isfield(evpOmega.parameters, char(fieldName)))
                testCase.verifyTrue(isfield(evpG.parameters, char(fieldName)))
                testCase.verifyTrue(isfield(evpF.parameters, char(fieldName)))
            end
            testCase.verifyTrue(isfield(evpK.parameters, "k"))
            testCase.verifyTrue(isfield(evpOmega.parameters, "omega"))
            testCase.verifyFalse(isfield(evpF.parameters, "k"))
            testCase.verifyFalse(isfield(evpK.parameters, "f0"))
            testCase.verifyFalse(isfield(evpK.parameters, "g"))
            testCase.verifyEqual(evpK.f0, f0, AbsTol=0)
            testCase.verifyEqual(evpK.g, g, AbsTol=0)

            lowLevelEVP = IMEigenvalueProblem(f0=f0, g=g, parameters=struct("f0", 3, "g", 4, "problemType", "custom"));
            testCase.verifyFalse(isfield(lowLevelEVP.parameters, "f0"))
            testCase.verifyFalse(isfield(lowLevelEVP.parameters, "g"))
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

        function basisSetDefaultsToPositionalModeIndexWhenOmitted(testCase)
            basisSet = IMBasisSet(nativeModes=zeros(3,4), h=1:4);

            testCase.verifyEqual(basisSet.modeIndex, 1:4)
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
            testCase.verifyEqual(basisG.modeIndex, 1:3)
            testCase.verifyEqual(basisF.modeIndex, 0:2)
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
                upperBoundary=IMBoundary.noSlip(), lowerBoundary=IMBoundary.rigid());
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
            positivePolicy = IMIndexPolicy.fromBoundaryRows(IMBoundaryRow.partialDepthPE(boundarySign="positive"));
            positiveIndex = positivePolicy.classify([1; 2; 3], struct());

            negativePolicy = IMIndexPolicy.fromBoundaryRows(IMBoundaryRow.partialDepthPE(boundarySign="negative"));
            negativeIndex = negativePolicy.classify([-2; -1; 3], struct());

            testCase.verifyEqual(positiveIndex.expectedNegativeCount, 0)
            testCase.verifyEqual(positiveIndex.negativeCount, 0)
            testCase.verifyEqual(negativeIndex.expectedNegativeCount, 2)
            testCase.verifyEqual(negativeIndex.negativeCount, 2)
        end

        function indexPolicyErrorsWhenObservedIndexDisagrees(testCase)
            policy = IMIndexPolicy.fromBoundaryRows(IMBoundaryRow.partialDepthPE(boundarySign="negative"));

            testCase.verifyError(@() policy.classify([-1; 2; 3], struct()), ...
                "IMIndexPolicy:IndexMismatch")
        end

        function indexPolicySelectsAndLabelsNamedModeBlocks(testCase)
            policy = IMIndexPolicy.fixed(expectedNegativeCount=2, expectedZeroCount=1, validationMode="none");

            selection = policy.selectModes([-10; -1; 0; 2; 5], 5, struct());

            testCase.verifyEqual(selection.sortIndex, [2 1 3 4 5])
            testCase.verifyEqual(selection.modeIndex, [-2 -1 0 1 2])
            testCase.verifyEqual(selection.index.negativeCount, 2)
            testCase.verifyEqual(selection.index.zeroCount, 1)
        end

        function indexPolicyExcludesUnexpectedNegativeArtifacts(testCase)
            policy = IMIndexPolicy.fixed(expectedNegativeCount=0, expectedZeroCount=1, validationMode="none");

            selection = policy.selectModes([-1e20; 0; 1; 2], 3, struct());

            testCase.verifyEqual(selection.sortIndex, [2 3 4])
            testCase.verifyEqual(selection.modeIndex, [0 1 2])
            testCase.verifyEqual(selection.index.negativeCount, 0)
            testCase.verifyEqual(selection.index.zeroCount, 1)
        end

        function hydrostaticFModeEVPDeclaresBarotropicZeroIndex(testCase)
            [~, ~, ~, ~, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);

            index = evp.indexPolicy.classify([0; 1; 2; 3], struct());

            testCase.verifyEqual(index.expectedZeroCount, 1)
            testCase.verifyEqual(index.zeroCount, 1)
            testCase.verifyEqual(index.positiveCount, 3)
        end

        function hydrostaticFModeEVPDoesNotDeclareWavenumber(testCase)
            [~, ~, ~, ~, g] = testCase.profile();
            evp = IMEigenvalueProblem.hydrostaticFModes(g=g);

            testCase.verifyFalse(isfield(evp.parameters, "k"))
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

            testCase.verifyTrue(isfield(evp.components.G, "operator"))
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

        function hydrostaticFModeGeostrophicNormalizationKeepsBarotropicNullMode(testCase)
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
            testCase.verifyEqual(basisSet.modeIndex, 0:4)
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
            testCase.verifyEqual(basisF.modeIndex, 0:4)
            testCase.verifyEqual(basisG.modeIndex, 1:4)
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
                upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.rigid());
            basisSet = solver.solveEVP(evp, nModes=3);
            z = solver.innerProductGrid(zDomain);
            rawG = solver.T*basisSet.nativeModes;
            weight = IMOperator.evaluateCoefficient(evp.components.G.innerWeight, z, testCase.evpContext(solver, evp));
            interior = solver.integrateInnerProduct(z, weight.*rawG(:,1).*rawG(:,1), zDomain);
            expectedFactor = sqrt(abs(interior + rawG(1,1)*rawG(1,1)));

            factors = basisSet.normalizationFactors(Normalization.unity);
            basisSet.normalization = Normalization.unity;
            fullGram = basisSet.gramMatrix("G");
            zCut = zDomain(2) - 1;
            partialGram = basisSet.partialGramMatrix("G", zDomain(1), zCut);
            GPartial = basisSet.evaluate("G", z);
            expectedPartial = solver.integrateInnerProduct(z, weight.*GPartial(:,1).*GPartial(:,1), [zDomain(1) zCut]);

            testCase.verifyEqual(factors(1), expectedFactor, RelTol=1e-10)
            testCase.verifyEqual(fullGram(1,1), 1, AbsTol=1e-10)
            testCase.verifyEqual(partialGram(1,1), expectedPartial, AbsTol=1e-10)
            testCase.verifyGreaterThan(fullGram(1,1), partialGram(1,1))
        end

        function fixedFrequencyKConstantUsesBoundaryInnerProductTerms(testCase)
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            zDomain = [-1300 0];
            nEVP = 48;
            f0 = 1e-4;
            g = 9.81;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtFrequency(omega=0.8*N0, f0=f0, g=g, ...
                upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.free());
            basisSet = solver.solveEVP(evp, nModes=3);
            z = solver.innerProductGrid(zDomain);
            rawG = solver.T*basisSet.nativeModes;
            weight = IMOperator.evaluateCoefficient(@(z,ctx) (ctx.N2(z) - ctx.f0*ctx.f0)/ctx.g, ...
                z, testCase.evpContext(solver, evp));
            interior = solver.integrateInnerProduct(z, weight.*rawG(:,1).*rawG(:,1), zDomain);
            expectedFactor = sqrt(abs(interior + rawG(1,1)*rawG(1,1) + rawG(end,1)*rawG(end,1)));

            factors = basisSet.normalizationFactors(Normalization.kConstant);

            testCase.verifyEqual(factors(1), expectedFactor, RelTol=1e-10)
        end

        function mixedInnerProductTermsContributeToGramMatrix(testCase)
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            zDomain = [-1300 0];
            nEVP = 48;
            g = 9.81;
            solver = IMSolverSpectral(N2=N2, zDomain=zDomain, nEVP=nEVP);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=0, g=g, ...
                upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.rigid());
            mixedTerms = [
                IMBoundaryRow.innerProductTerm("G", "surface", 1, IMBoundaryRow.trace("G"), IMBoundaryRow.trace("F"))
                IMBoundaryRow.innerProductTerm("G", "surface", 1, IMBoundaryRow.trace("F"), IMBoundaryRow.trace("G"))
            ];
            evp.boundaryRows(end+1,1) = IMBoundaryRow.active(endpoint="upper", component="G", ...
                indexSign=1, innerProductTerms=mixedTerms);
            basisSet = solver.solveEVP(evp, nModes=3);
            basisSet.normalization = Normalization.unity;

            z = solver.innerProductGrid(zDomain);
            weight = IMOperator.evaluateCoefficient(evp.components.G.innerWeight, z, testCase.evpContext(solver, evp));
            G = basisSet.evaluate("G", z);
            interior = solver.integrateInnerProduct(z, weight.*G(:,1).*G(:,1), zDomain);
            GSurface = basisSet.evaluate("G", zDomain(2));
            FSurface = basisSet.evaluate("F", zDomain(2));
            expected = interior + GSurface(1,1)*GSurface(1,1) ...
                + GSurface(1,1)*FSurface(1,1) + FSurface(1,1)*GSurface(1,1);

            gram = basisSet.gramMatrix("G");

            testCase.verifyEqual(gram(1,1), expected, AbsTol=1e-10)
        end

        function linearBoundaryFamiliesDeclareTrustedInnerProductTerms(testCase)
            linearF = IMBoundary.linearF(c=2, d=3).resolve(endpoint="upper", primaryComponent="F");
            linearG = IMBoundary.linearG(a=2, b=3).resolve(endpoint="upper", primaryComponent="G");

            testCase.verifyTrue(linearF.hasKnownInnerProduct)
            testCase.verifyEqual(linearF.innerProductTerms(1).innerProductComponent, "G")
            testCase.verifyEqual(linearF.innerProductTerms(1).location, "surface")
            testCase.verifyEqual(linearF.innerProductTerms(1).coefficient, -1.5, AbsTol=0)
            testCase.verifyTrue(linearG.hasKnownInnerProduct)
            testCase.verifyEqual(linearG.innerProductTerms(1).innerProductComponent, "G")
            testCase.verifyEqual(linearG.innerProductTerms(1).coefficient, 1.5, AbsTol=0)
        end

        function unresolvedLinearFamiliesWarnButStillAssembleRows(testCase)
            testCase.verifyWarning(@() IMBoundary.linearF(a=1, b=2, c=3).resolve(endpoint="upper", primaryComponent="F"), ...
                "IMBoundary:UnknownInnerProduct")
            warningState = warning("off", "IMBoundary:UnknownInnerProduct");
            cleanup = onCleanup(@() warning(warningState));
            boundary = IMBoundary.linearF(a=1, b=2, c=3).resolve(endpoint="upper", primaryComponent="F");

            testCase.verifyFalse(boundary.hasKnownInnerProduct)
            testCase.verifyEmpty(boundary.innerProductTerms)
            testCase.verifyError(@() IMBoundary.linearF(c=1).resolve(endpoint="upper", primaryComponent="G"), ...
                "IMBoundary:UnsupportedResolution")
        end

        function finiteDifferenceUnityNormalizationUsesBoundedTrapz(testCase)
            [N2, zDomain, ~, f0, g] = testCase.profile();
            solver = IMSolverFiniteDifference(z=linspace(zDomain(1), zDomain(2), 35).', N2=N2);
            evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=f0, g=g, ...
                upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.rigid());
            basisSet = solver.solveEVP(evp, nModes=3);
            z = solver.innerProductGrid(zDomain);
            rawG = solver.evaluateNativeModes(basisSet.nativeModes, z);
            weight = IMOperator.evaluateCoefficient(evp.components.G.innerWeight, z, testCase.evpContext(solver, evp));
            normValue = trapz(z, weight.*rawG(:,1).*rawG(:,1));
            surfaceValue = solver.evaluateNativeModes(basisSet.nativeModes(:,1), zDomain(2));
            expectedFactor = sqrt(abs(normValue + surfaceValue*surfaceValue));

            factors = basisSet.normalizationFactors(Normalization.unity);

            testCase.verifyEqual(factors(1), expectedFactor, RelTol=1e-10)
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
            testCase.verifyEqual(basisSet.modeIndex, 1:nModes)
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
                upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.rigid());
            basisSet = IMBasisSet.constantStratification(evp=evp, N0=N0, zDomain=zDomain, nModes=nModes);

            direct = InternalModesConstantStratification(N0=N0, zIn=zDomain, zOut=z, latitude=0, nModes=nModes, g=g);
            direct.upperBoundary = UpperBoundary.freeSurface;
            direct.normalization = Normalization.kConstant;
            [FExpected, GExpected, hExpected] = direct.modesAtWavenumber(k);

            testCase.verifyEqual(basisSet.h, hExpected, RelTol=1e-10)
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
            testCase.verifyError(@() IMBasisSet.exponentialStratification(), ...
                "IMBasisSet:AnalyticalBasisNotImplemented")
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
            weight = IMOperator.evaluateCoefficient(evp.components.G.innerWeight, z, testCase.evpContext(solver, evp));
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

        function [surfaceWeight, bottomWeight] = scalarEndpointWeights(~, evp, component)
            surfaceWeight = 0;
            bottomWeight = 0;
            for iBoundary = 1:length(evp.boundaryRows)
                terms = evp.boundaryRows(iBoundary).innerProductTerms;
                for iTerm = 1:length(terms)
                    term = terms(iTerm);
                    if string(term.innerProductComponent) ~= string(component)
                        continue;
                    end
                    if ~isnumeric(term.coefficient) || ~isscalar(term.coefficient)
                        continue;
                    end
                    if string(term.leftTrace.component) ~= string(component) ...
                            || string(term.rightTrace.component) ~= string(component)
                        continue;
                    end
                    if term.leftTrace.derivativeOrder ~= 0 || term.rightTrace.derivativeOrder ~= 0
                        continue;
                    end
                    if string(term.location) == "surface"
                        surfaceWeight = surfaceWeight + term.coefficient;
                    elseif string(term.location) == "bottom"
                        bottomWeight = bottomWeight + term.coefficient;
                    end
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
