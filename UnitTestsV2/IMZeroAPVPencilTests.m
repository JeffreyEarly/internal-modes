classdef IMZeroAPVPencilTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function repeatedIndefinitePencilHasACompleteSignedOrthonormalFrame(testCase)
            right = [2 1;1 -2];
            [C,lambda,signatures,evidence] = solve(right,right);
            testCase.verifyEqual(lambda,[1;1],AbsTol=3e-14);
            testCase.verifyEqual(C'*right*C,diag(signatures),AbsTol=3e-14);
            testCase.verifyEqual(sort(signatures),[-1;1]);
            testCase.verifyLessThan(evidence.rightOrthogonalityError,3e-14);
            testCase.verifyLessThan(max(evidence.relativePencilResidual),3e-14);
        end

        function nullCoordinateEigenvectorsDoNotInvalidateANondegenerateRepeatedSpace(testCase)
            right = [0 1;1 0];
            [C,lambda,signatures] = solve(right,right);
            testCase.verifyEqual(lambda,[1;1],AbsTol=3e-14);
            testCase.verifyEqual(C'*right*C,diag(signatures),AbsTol=3e-14);
            testCase.verifyEqual(sort(signatures),[-1;1]);
            testCase.verifyLessThan(cond(C),1+3e-14);
        end

        function uniformPhysicalUnitsDoNotChangeMetricRankOrEigenvalues(testCase)
            left = [3 .2;.2 -1]; right = [2 .3;.3 -2];
            [~,reference] = solve(left,right);
            for scale = [1e-250 1 1e250]
                [C,lambda,signatures,evidence] = solve(scale*left,scale*right);
                testCase.verifyEqual(lambda,reference,RelTol=3e-13,AbsTol=3e-14);
                testCase.verifyEqual(C'*(scale*right)*C,diag(signatures),AbsTol=3e-13);
                testCase.verifyLessThan(max(evidence.relativePencilResidual),3e-13);
            end
        end

        function zeroLeftMatrixRetainsTheEntireSignedSpace(testCase)
            right = 1e-250*[0 1;1 0];
            [C,lambda,signatures,evidence] = solve(zeros(2),right);
            testCase.verifyEqual(lambda,zeros(2,1),AbsTol=0);
            testCase.verifyEqual(C'*right*C,diag(signatures),AbsTol=3e-14);
            testCase.verifyEqual(evidence.relativePencilResidual,zeros(2,1),AbsTol=0);
        end

        function smallComplexSpectrumCannotBeCoercedToReal(testCase)
            testCase.verifyError(@() solve(1e-220*[0 1;1 0],diag([1 -1])),"IMZeroAPVPencilAudit:InvalidPencilSpectrum");
        end

        function singularAndDefectivePencilsFailExplicitly(testCase)
            testCase.verifyError(@() solve(eye(2),diag([1 0])),"IMZeroAPVPencilAudit:NullRightMetric");
            testCase.verifyError(@() solve(eye(2),zeros(2)),"IMZeroAPVPencilAudit:NullRightMetric");
            testCase.verifyError(@() solve(ones(2),diag([1 -1])),"IMZeroAPVPencilAudit:SingularCoordinateRotation");
        end

        function scaleRestorationAvoidsOverflowButRejectsUnrepresentableSpectra(testCase)
            left = 1e308*eye(2); right = 1e308*diag([.1 -1]);
            [C,lambda,signatures,evidence] = solve(left,right);
            testCase.verifyEqual(lambda,[10;-1],RelTol=3e-14);
            testCase.verifyEqual(C'*right*C,diag(signatures),AbsTol=3e-14);
            testCase.verifyLessThan(max(evidence.relativePencilResidual),3e-14);
            testCase.verifyError(@() solve(1e308*eye(2),1e-308*eye(2)),"IMZeroAPVPencilAudit:InvalidPencilSpectrum");
            testCase.verifyError(@() solve(1e-308*eye(2),1e308*eye(2)),"IMZeroAPVPencilAudit:InvalidPencilSpectrum");
            testCase.verifyError(@() solve(diag([1e308 1e-100]),1e308*eye(2)),"IMZeroAPVPencilAudit:InvalidPencilSpectrum");
        end

        function nearlyRepeatedRootsStillSatisfyTheirPencilAndMetric(testCase)
            right = [2 1;1 -2];
            [directions,metric] = eig(right,"vector");
            for separation = [1e-14 1e-7]
                left = directions*diag(metric.*[1;1+separation])*directions';
                [C,lambda,signatures,evidence] = solve(left,right);
                testCase.verifyEqual(sort(lambda),[1;1+separation],RelTol=3e-13);
                testCase.verifyEqual(C'*right*C,diag(signatures),AbsTol=3e-13);
                testCase.verifyLessThan(max(evidence.relativePencilResidual),1e-10);
            end
        end

        function emptyEndpointPagesHaveAnEmptyIdentityWithoutAnEigensolve(testCase)
            [C,lambda,signatures,evidence] = IMGeostrophicZeroAPVFormTools.solveSymmetricPencil(zeros(0,0,2),zeros(0,0,2),0,2,"descending","IMZeroAPVPencilAudit");
            testCase.verifySize(C,[0 0 2]);
            testCase.verifySize(lambda,[0 2]);
            testCase.verifySize(signatures,[0 2]);
            testCase.verifyEqual(evidence.rightOrthogonalityError,[0 0],AbsTol=0);
            testCase.verifyEqual(evidence.rotationConditionNumber,[1 1],AbsTol=0);
            IMGeostrophicZeroAPVFormTools.validateRotationPages(C,0,2,"IMZeroAPVPencilAudit");
        end

        function shapeAndSymmetryValidationUseRelativeMatrixUnits(testCase)
            prefix = "IMZeroAPVPencilAudit";
            testCase.verifyError(@() IMGeostrophicZeroAPVFormTools.expandMatrixPages(ones(2,2,1,2),2,1,"matrix",prefix),prefix+":InvalidPencilShape");
            testCase.verifyError(@() IMGeostrophicZeroAPVFormTools.expandMatrixPages(1e-250*[1 1;0 1],2,1,"matrix",prefix),prefix+":NonSymmetricPencil");
            testCase.verifyError(@() IMGeostrophicZeroAPVFormTools.expandMatrixPages(1i*eye(2),2,1,"matrix",prefix),prefix+":InvalidPencilMatrix");
            testCase.verifyError(@() IMGeostrophicZeroAPVFormTools.expandMatrixPages(diag([1 Inf]),2,1,"matrix",prefix),prefix+":InvalidPencilMatrix");
            testCase.verifyFalse(IMGeostrophicZeroAPVFormTools.pagesEqual(1e-250*eye(2),2e-250*eye(2)));
            testCase.verifyEqual(IMGeostrophicZeroAPVFormTools.symmetrize(1e308*eye(2)),1e308*eye(2),AbsTol=0);
        end
    end
end

function [C,lambda,signatures,evidence] = solve(left,right)
[C,lambda,signatures,evidence] = IMGeostrophicZeroAPVFormTools.solveSymmetricPencil(left,right,2,1,"descending","IMZeroAPVPencilAudit");
end
