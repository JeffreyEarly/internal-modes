classdef IMQuadraturePointTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
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
        function dirichletGridUsesNextSelectedMode(testCase)
            [solver, evp] = testCase.regularProblem(96);
            basisWithoutAuxiliary = solver.solveEVP(evp,nModes=4);
            basisWithAuxiliary = solver.solveEVP(evp,nModes=5);

            zAutomatic = basisWithoutAuxiliary.quadraturePoints();
            zRetained = basisWithAuxiliary.quadraturePoints(nModes=4);
            expected = linspace(-1,0,6).';

            testCase.verifyEqual(zAutomatic,expected,AbsTol=1e-10)
            testCase.verifyEqual(zRetained,zAutomatic,AbsTol=1e-12)
            testCase.verifyEqual(zAutomatic(1),-1,AbsTol=0)
            testCase.verifyEqual(zAutomatic(end),0,AbsTol=0)
        end

        function retainedPrefixIncludesNullMode(testCase)
            solver = IMSolverSpectral(nEVP=96);
            N2 = @(z) ones(size(z));
            evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-1 0]);
            basisSet = solver.solveEVP(evp,nModes=5);
            z = basisSet.quadraturePoints(nModes=4);
            F = basisSet.F(z(2:end-1));

            testCase.verifyEqual(basisSet.modeNumber,0:4,AbsTol=0)
            testCase.verifyEqual(length(z),6)
            testCase.verifyLessThan(max(abs(F(:,5))),1e-11)
        end

        function activeRobinGridKeepsPhysicalEndpoints(testCase)
            solver = IMSolverSpectral(nEVP=96);
            surfaceBoundary = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1, ...
                surfaceBoundary=surfaceBoundary,bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=5);
            z = basisSet.quadraturePoints(nModes=4);
            values = basisSet.u(z);

            testCase.verifyEqual(z(1),-1,AbsTol=0)
            testCase.verifyEqual(z(end),0,AbsTol=0)
            testCase.verifyGreaterThan(abs(values(end,5)),1e-3)
            testCase.verifyLessThan(max(abs(values(2:end-1,5))),1e-10)
        end

        function allSpectralCoordinatesRecoverPhysicalRoots(testCase)
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-4000 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            coordinateKinds = ["z" "wkb" "density"];
            grids = zeros(6,length(coordinateKinds));
            for iKind = 1:length(coordinateKinds)
                solver = IMSolverSpectral(nEVP=128,coordinateKind=coordinateKinds(iKind));
                evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain);
                basisSet = solver.solveEVP(evp,nModes=5);
                grids(:,iKind) = basisSet.quadraturePoints(nModes=4);
                G = basisSet.G(grids(2:end-1,iKind));
                testCase.verifyLessThan(max(abs(G(:,5))),1e-8)
            end

            testCase.verifyEqual(grids(:,2),grids(:,1),AbsTol=1e-3)
            testCase.verifyEqual(grids(:,3),grids(:,1),AbsTol=1e-3)
        end

        function generatedGridBuildsFittedTransform(testCase)
            [solver, evp] = testCase.regularProblem(96);
            basisSet = solver.solveEVP(evp,nModes=4);
            z = basisSet.quadraturePoints();
            [weights, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=4);
            transform = basisSet.discreteTransform(z=z,nModes=4);

            testCase.verifyGreaterThan(fit.exitFlag,0)
            testCase.verifyEqual(weights,fit.weights,AbsTol=0)
            testCase.verifyGreaterThanOrEqual(min(fit.weights),-1e-12)
            testCase.verifyEqual(sum(fit.weights),1,AbsTol=1e-10)
            testCase.verifyEqual(transform.z,z,AbsTol=0)
            testCase.verifyLessThan(transform.roundTripError,1e-11)
        end

        function unsupportedAndInvalidRequestsThrowStructuredErrors(testCase)
            [solver, evp] = testCase.regularProblem(96);
            basisSet = solver.solveEVP(evp,nModes=4);
            testCase.verifyError(@() basisSet.quadraturePoints(nModes=5),"IMBasisSet:InvalidQuadratureModeCount")

            lowResolutionSolver = IMSolverSpectral(nEVP=4);
            lowResolutionBasis = lowResolutionSolver.solveEVP(evp,nModes=2);
            testCase.verifyError(@() lowResolutionBasis.quadraturePoints(),"IMBasisSet:AuxiliaryModeUnavailable")

            inconsistentBasis = IMBasisSet(solver=basisSet.solver,evp=evp,nativeModes=basisSet.nativeModes, ...
                eigenvalues=basisSet.eigenvalues+1,modeNumber=basisSet.modeNumber);
            testCase.verifyError(@() inconsistentBasis.quadraturePoints(),"IMBasisSet:AuxiliaryModeMismatch")

            finiteDifferenceSolver = IMSolverFiniteDifference(z=linspace(-1,0,65).');
            finiteDifferenceBasis = finiteDifferenceSolver.solveEVP(evp,nModes=3);
            testCase.verifyError(@() finiteDifferenceBasis.quadraturePoints(nModes=2),"IMSolver:UnsupportedModeRoots")
        end

        function deficientRootSetIsRejected(testCase)
            [solver, evp] = testCase.regularProblem(8);
            solver = solver.configuredForEVP(evp);
            nativeModes = zeros(solver.nEVP,4);
            nativeModes(1,:) = 1;
            basisSet = IMBasisSet(solver=solver,evp=evp,nativeModes=nativeModes, ...
                eigenvalues=1:4,modeNumber=1:4);

            testCase.verifyError(@() basisSet.quadraturePoints(nModes=3),"IMBasisSet:InsufficientQuadraturePoints")
        end
    end

    methods (Access = private)
        function [solver, evp] = regularProblem(~,nEVP)
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1, ...
                surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
        end
    end
end
