classdef IMDiscreteTransformAssessmentTests < matlab.unittest.TestCase

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
        function exactPointConstructionCoversStratificationsAndCoordinates(testCase)
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            coordinateKinds = ["z" "wkb" "density"];
            constantN2 = @(z) N0*N0*ones(size(z));
            exponentialN2 = @(z) N0*N0*exp(2*z/1300);
            evps = [IMInternalModes.hydrostaticFModes(N2=constantN2,zDomain=zDomain), ...
                IMInternalModes.hydrostaticFModes(N2=exponentialN2,zDomain=zDomain), ...
                IMInternalModes.hydrostaticGModes(N2=exponentialN2,zDomain=zDomain)];

            for iEVP = 1:length(evps)
                for iCoordinate = 1:length(coordinateKinds)
                    solver = IMSolverSpectral(nEVP=128,coordinateKind=coordinateKinds(iCoordinate));
                    basisSet = solver.solveEVP(evps(iEVP),nModes=12);
                    [transform,assessment] = basisSet.discreteTransform(nPoints=8);

                    testCase.verifyClass(assessment,"IMDiscreteTransformAssessment")
                    testCase.verifyEqual(assessment.requestedPointCount,8)
                    testCase.verifyEqual(assessment.actualPointCount,8)
                    testCase.verifyEqual(length(transform.z),8)
                    testCase.verifyEqual(transform.z(1),zDomain(1),AbsTol=0)
                    testCase.verifyEqual(transform.z(end),zDomain(2),AbsTol=0)
                    testCase.verifyEqual(assessment.candidateModeCount,6)
                    testCase.verifyEqual(assessment.retainedModeCount,6)
                    testCase.verifyEqual(assessment.gramPolicy.tolerance,1e-2,AbsTol=0)
                    testCase.verifyLessThanOrEqual(assessment.prefixDiagnostics.gramError(end),1e-2)
                end
            end
        end

        function assessmentReproducesFixedRuleAndNestedFit(testCase)
            basisSet = testCase.constantFBasis(10);
            z = -4000 + 4000*linspace(0,1,17).'.^1.25;
            [transform,assessment] = basisSet.discreteTransform(z=z,nModes=6);
            candidate = assessment.candidateTransform;
            fit = assessment.weightFit;
            finalRow = assessment.prefixDiagnostics(end,:);

            testCase.verifyClass(fit,"IMQuadratureWeightFit")
            testCase.verifyEqual(fit.transform.weights,candidate.weights,AbsTol=0)
            testCase.verifyEqual(fit.transform.forwardMatrix,candidate.forwardMatrix,AbsTol=0)
            testCase.verifyEqual(fit.transform.inverseMatrix,candidate.inverseMatrix,AbsTol=0)
            testCase.verifyEqual(transform.weights,candidate.weights,AbsTol=0)
            testCase.verifyEqual(assessment.transform.forwardMatrix,transform.forwardMatrix,AbsTol=0)
            testCase.verifyEqual(assessment.candidateModeNumber,basisSet.modeNumber(1:6),AbsTol=0)
            testCase.verifyEqual(assessment.retainedModeNumber,transform.modeNumber,AbsTol=0)
            testCase.verifyEqual(finalRow.gramError,candidate.relativeGramOperatorError,AbsTol=0)
            testCase.verifyEqual(finalRow.roundTripError,candidate.roundTripError,AbsTol=0)
            testCase.verifyEqual(finalRow.inverseMatrixConditionNumber,candidate.inverseMatrixConditionNumber,AbsTol=0)
            testCase.verifyEqual(finalRow.gramConditionNumber,candidate.gramConditionNumber,AbsTol=0)
            testCase.verifyEqual(height(assessment.prefixDiagnostics),assessment.candidateModeCount)
            testCase.verifyEqual(assessment.limitingPolicy,"none")
        end

        function suppliedWeightsBypassOnlyTheFit(testCase)
            basisSet = testCase.constantFBasis(8);
            z = linspace(-4000,0,17).';
            weights = testCase.trapezoidalWeights(z);
            [transform,assessment] = basisSet.discreteTransform(z=z,weights=weights,nModes=6);

            testCase.verifyEmpty(assessment.weightFit)
            testCase.verifyEqual(transform.weights,weights,AbsTol=0)
            testCase.verifyEqual(assessment.candidateTransform.weights,weights,AbsTol=0)
            testCase.verifyEqual(assessment.requestedPointCount,length(z))
            testCase.verifyEqual(assessment.actualPointCount,length(z))
        end

        function unattainablePointCountsReportNearbyCounts(testCase)
            basisSet = testCase.constantFBasis(12);
            try
                basisSet.discreteTransform(nPoints=15);
                testCase.assertFail("Expected an unattainable exact point-count error.")
            catch exception
                testCase.verifyEqual(string(exception.identifier),"IMBasisSet:UnattainableDiscretePointCount")
                testCase.verifySubstring(exception.message,"exactly 15 points")
                testCase.verifySubstring(exception.message,"14")
            end
        end

        function pointAndBandSpecificationsAreUnambiguous(testCase)
            basisSet = testCase.constantFBasis(8);
            z = linspace(-4000,0,9).';
            weights = testCase.trapezoidalWeights(z);

            testCase.verifyError(@() basisSet.discreteTransform(),"IMBasisSet:InvalidDiscretePointSpecification")
            testCase.verifyError(@() basisSet.discreteTransform(nPoints=8,z=z),"IMBasisSet:InvalidDiscretePointSpecification")
            testCase.verifyError(@() basisSet.discreteTransform(weights=weights),"IMBasisSet:InvalidDiscretePointSpecification")
            testCase.verifyError(@() basisSet.discreteTransform(nPoints=8,nModes=6),"IMBasisSet:InvalidDiscretePointSpecification")
            testCase.verifyError(@() basisSet.discreteTransform(nPoints=[7 8]),"IMBasisSet:InvalidDiscreteTransformOption")
        end
    end

    methods (Access = private)
        function basisSet = constantFBasis(~,nModes)
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            N2 = @(z) N0*N0*ones(size(z));
            solver = IMSolverSpectral(nEVP=128);
            evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=zDomain);
            basisSet = solver.solveEVP(evp,nModes=nModes);
        end

        function weights = trapezoidalWeights(~,z)
            weights = [(z(2)-z(1))/2; (z(3:end)-z(1:end-2))/2; (z(end)-z(end-1))/2];
        end
    end
end
