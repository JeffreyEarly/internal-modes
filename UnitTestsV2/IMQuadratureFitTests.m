classdef IMQuadratureFitTests < matlab.unittest.TestCase

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
        function omittedIncrementsUseDefaultPhysicalFit(testCase)
            basisSet = testCase.regularBasis(4);
            z = -1 + linspace(0,1,17).'.^1.4;
            fit = basisSet.fitQuadrature(z=z,nModes=3);
            transform = basisSet.discreteTransform(z=z,nModes=3);

            testCase.verifyClass(fit,"IMQuadratureFit")
            testCase.verifyEqual(transform.increments,fit.fittedIncrements,AbsTol=1e-12)
            testCase.verifyEqual(transform.forwardMatrix,fit.fittedTransform.forwardMatrix,RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(fit.objectiveName,"normalizedGram")
            testCase.verifyTrue(fit.nonnegativeConstraint)
            testCase.verifyTrue(fit.depthConstraint)
            testCase.verifyGreaterThan(fit.exitFlag,0)
            testCase.verifyGreaterThanOrEqual(min(fit.fittedIncrements),-1e-12)
            testCase.verifyEqual(sum(fit.fittedIncrements),1,AbsTol=1e-10)
            testCase.verifyLessThanOrEqual(fit.fittedResidualNorm,fit.geometricResidualNorm + 1e-10)
        end

        function suppliedIncrementsPreservePhaseOnePath(testCase)
            basisSet = testCase.regularBasis(3);
            z = linspace(-1,0,21).';
            dz = testCase.geometricIncrements(z,[-1 0]);
            transform = basisSet.discreteTransform(z=z,increments=dz,nModes=2);

            testCase.verifyEqual(transform.increments,dz,AbsTol=0)
            testCase.verifyEqual(transform.metricMatrix,diag(dz),AbsTol=0)
        end

        function geometricIncrementsUseFullDomainControlVolumes(testCase)
            basisSet = testCase.regularBasis(3);
            zInterior = [-0.9; -0.63; -0.31; -0.08];
            fitInterior = basisSet.fitQuadrature(z=zInterior,nModes=2);
            expectedInterior = testCase.geometricIncrements(zInterior,[-1 0]);

            zEndpoints = [-1; -0.72; -0.4; -0.15; 0];
            fitEndpoints = basisSet.fitQuadrature(z=zEndpoints,nModes=2);
            expectedEndpoints = testCase.geometricIncrements(zEndpoints,[-1 0]);

            testCase.verifyEqual(fitInterior.geometricIncrements,expectedInterior,AbsTol=1e-14)
            testCase.verifyEqual(fitEndpoints.geometricIncrements,expectedEndpoints,AbsTol=1e-14)
            testCase.verifyEqual(sum(fitInterior.geometricIncrements),1,AbsTol=1e-14)
            testCase.verifyEqual(sum(fitEndpoints.geometricIncrements),1,AbsTol=1e-14)
        end

        function normalizedGramSystemIncludesFixedEndpointContribution(testCase)
            zDomain = [-1 0];
            solver = IMSolverSpectral(nEVP=96);
            surfaceBoundary = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            evp = IMEigenvalueProblem(zDomain=zDomain,p=1,q=0,r=1,surfaceBoundary=surfaceBoundary,bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=3);
            z = linspace(-1,0,17).';
            fit = basisSet.fitQuadrature(z=z,nModes=2);

            Phi = fit.geometricTransform.basisMatrix;
            target = fit.geometricTransform.targetGramMatrix;
            targetNorms = diag(target);
            endpointWeight = evp.innerProduct().surfaceWeights(1);
            endpointGram = endpointWeight.coefficient*endpointWeight.c^2*(Phi(end,:).'*Phi(end,:));
            expectedA = zeros(4,length(z));
            expectedB = zeros(4,1);
            iRow = 0;
            for iMode = 1:2
                for jMode = 1:2
                    iRow = iRow + 1;
                    scale = sqrt(abs(targetNorms(iMode)*targetNorms(jMode)));
                    expectedA(iRow,:) = (Phi(:,iMode).*Phi(:,jMode)).'/scale;
                    expectedB(iRow) = (target(iMode,jMode) - endpointGram(iMode,jMode))/scale;
                end
            end

            testCase.verifyEqual(fit.objectiveMatrix,expectedA,RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(fit.objectiveTarget,expectedB,RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(fit.fittedObjectiveResidual,expectedA*fit.fittedIncrements - expectedB,AbsTol=1e-13)
            testCase.verifyEqual(fit.geometricObjectiveResidual,expectedA*fit.geometricIncrements - expectedB,AbsTol=1e-13)
            testCase.verifyError(@() basisSet.fitQuadrature(z=z(1:end-1),nModes=2),"IMBasisSet:MissingDiscreteEndpointSample")
        end

        function derivativeEndpointMetricRemainsUnsupported(testCase)
            solver = IMSolverSpectral(nEVP=96);
            surfaceBoundary = IMBoundaryCondition(a=1,b=0,c=0,d=1);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=surfaceBoundary,bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=2);
            z = linspace(-1,0,17).';

            testCase.verifyError(@() basisSet.fitQuadrature(z=z),"IMBasisSet:UnsupportedDiscreteEndpointMetric")
        end

        function fittedRuleImprovesVariableWeightGramDiagnostic(testCase)
            zDomain = [-1 0];
            solver = IMSolverSpectral(nEVP=128);
            evp = IMEigenvalueProblem(zDomain=zDomain,p=1,q=0,r=@(z) exp(2*z),surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=6);
            z = -1 + linspace(0,1,15).'.^1.7;
            fit = basisSet.fitQuadrature(z=z,nModes=5);

            testCase.verifyLessThan(fit.fittedTransform.relativeGramError,fit.geometricTransform.relativeGramError)
            testCase.verifyLessThanOrEqual(fit.fittedResidualNorm,fit.geometricResidualNorm + 1e-10)
        end

        function customObjectiveCanReplaceNormalizedGramSystem(testCase)
            basisSet = testCase.regularBasis(3);
            z = [-1; -0.73; -0.44; -0.18; 0];
            desiredIncrements = [0.08; 0.17; 0.24; 0.29; 0.22];
            objective = @(context) struct("A",eye(length(context.z)),"b",desiredIncrements,"name","prescribedIncrements");
            fit = basisSet.fitQuadrature(z=z,nModes=2,objective=objective);

            testCase.verifyEqual(fit.objectiveName,"prescribedIncrements")
            testCase.verifyEqual(fit.objectiveMatrix,eye(length(z)),AbsTol=0)
            testCase.verifyEqual(fit.objectiveTarget,desiredIncrements,AbsTol=0)
            testCase.verifyEqual(fit.fittedIncrements,desiredIncrements,AbsTol=1e-10)
        end

        function customObjectiveCanReweightAndAppendRows(testCase)
            basisSet = testCase.regularBasis(3);
            z = -1 + linspace(0,1,9).'.^1.3;
            baseline = basisSet.fitQuadrature(z=z,nModes=2);
            expectedA = [2*baseline.objectiveMatrix; eye(length(z))];
            expectedB = [2*baseline.objectiveTarget; baseline.geometricIncrements];
            objective = @(context) struct("A",[2*context.normalizedGramA; eye(length(context.z))], ...
                "b",[2*context.normalizedGramB; context.geometricIncrements],"name","weightedRegularizedGram");
            fit = basisSet.fitQuadrature(z=z,nModes=2,objective=objective);

            testCase.verifyEqual(fit.objectiveName,"weightedRegularizedGram")
            testCase.verifyEqual(fit.objectiveMatrix,expectedA,RelTol=1e-13,AbsTol=1e-13)
            testCase.verifyEqual(fit.objectiveTarget,expectedB,RelTol=1e-13,AbsTol=1e-13)
        end

        function constraintsCanBeDisabledIndependently(testCase)
            basisSet = testCase.regularBasis(3);
            z = [-1; -0.75; -0.5; -0.25; 0];
            signedIncrements = [-0.02; 0.2; 0.3; 0.3; 0.22];
            signedObjective = @(context) struct("A",eye(length(context.z)),"b",signedIncrements,"name","signed");
            signedFit = basisSet.fitQuadrature(z=z,nModes=2,objective=signedObjective,nonnegative=false);

            shallowIncrements = 0.75*testCase.geometricIncrements(z,[-1 0]);
            shallowObjective = @(context) struct("A",eye(length(context.z)),"b",shallowIncrements,"name","shallow");
            shallowFit = basisSet.fitQuadrature(z=z,nModes=2,objective=shallowObjective,constrainDepth=false);

            testCase.verifyFalse(signedFit.nonnegativeConstraint)
            testCase.verifyTrue(signedFit.depthConstraint)
            testCase.verifyEqual(signedFit.fittedIncrements,signedIncrements,AbsTol=1e-10)
            testCase.verifyTrue(signedFit.fittedTransform.hasNegativeIncrements)
            testCase.verifyEqual(sum(signedFit.fittedIncrements),1,AbsTol=1e-10)
            testCase.verifyTrue(shallowFit.nonnegativeConstraint)
            testCase.verifyFalse(shallowFit.depthConstraint)
            testCase.verifyEqual(shallowFit.fittedIncrements,shallowIncrements,AbsTol=1e-10)
            testCase.verifyEqual(sum(shallowFit.fittedIncrements),0.75,AbsTol=1e-10)
        end

        function malformedObjectivesThrowStructuredErrors(testCase)
            basisSet = testCase.regularBasis(3);
            z = linspace(-1,0,7).';
            missingTarget = @(context) struct("A",eye(length(context.z)));
            wrongColumns = @(context) struct("A",ones(2,length(context.z) + 1),"b",ones(2,1));
            nonfinite = @(context) struct("A",NaN(1,length(context.z)),"b",1);
            failing = @(context) error("Test:ObjectiveFailure","failed for %d points",length(context.z));

            testCase.verifyError(@() basisSet.fitQuadrature(z=z,nModes=2,objective="notAnObjective"),"IMBasisSet:UnknownQuadratureObjective")
            testCase.verifyError(@() basisSet.fitQuadrature(z=z,nModes=2,objective=missingTarget),"IMBasisSet:InvalidQuadratureObjective")
            testCase.verifyError(@() basisSet.fitQuadrature(z=z,nModes=2,objective=wrongColumns),"IMBasisSet:InvalidQuadratureObjective")
            testCase.verifyError(@() basisSet.fitQuadrature(z=z,nModes=2,objective=nonfinite),"IMBasisSet:InvalidQuadratureObjective")
            testCase.verifyError(@() basisSet.fitQuadrature(z=z,nModes=2,objective=failing),"IMBasisSet:QuadratureObjectiveFailed")
        end
    end

    methods (Access = private)
        function basisSet = regularBasis(~,nModes)
            solver = IMSolverSpectral(nEVP=96);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=nModes);
        end

        function dz = geometricIncrements(~,z,zDomain)
            edges = [zDomain(1); 0.5*(z(1:end-1) + z(2:end)); zDomain(2)];
            dz = diff(edges);
        end
    end
end
