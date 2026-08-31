classdef IMSolverIntegrationWeightsTests < matlab.unittest.TestCase

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
        function restorePath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function fullDomainSpectralWeightsMatchScalarIntegrator(testCase)
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-D 0]);
            basis = IMSolverSpectral(nEVP=128,coordinateKind="wkb").solveEVP(evp,nModes=4);
            solver = basis.solver;
            z = solver.innerProductGrid(evp.zDomain);
            values = [sin(z/317) cos(z/541) 1+(z/D).^3];
            weights = solver.innerProductWeights(z,evp.zDomain);

            for iValue = 1:size(values,2)
                expected = solver.integrateInnerProduct(z,values(:,iValue),evp.zDomain);
                testCase.verifyEqual(weights.'*values(:,iValue),expected,RelTol=5e-12,AbsTol=5e-12)
            end
        end

        function partialSpectralWeightsUseExactBoundedOperator(testCase)
            zDomain = [-4000 0];
            partialBounds = [-3150 -275];
            evp = IMEigenvalueProblem(zDomain=zDomain,p=1,q=0,r=1);
            basis = IMSolverSpectral(nEVP=96).solveEVP(evp,nModes=4);
            solver = basis.solver;
            z = solver.innerProductGrid(partialBounds);
            values = exp(z/2300).*cos(z/411);
            weights = solver.innerProductWeights(z,partialBounds);
            expected = solver.integrateInnerProduct(z,values,partialBounds);

            testCase.verifyEqual(weights.'*values,expected,RelTol=2e-13,AbsTol=2e-13)
        end

        function finiteDifferenceWeightsMatchTrapezoidalIntegrator(testCase)
            zDomain = [-4000 0];
            zNative = -4000+4000*linspace(0,1,83).'.^1.3;
            evp = IMEigenvalueProblem(zDomain=zDomain,p=1,q=0,r=1);
            solver = IMSolverFiniteDifference(z=zNative).configuredForEVP(evp);
            z = solver.innerProductGrid(zDomain);
            values = sin(z/317)+0.2*cos(z/901);
            weights = solver.innerProductWeights(z,zDomain);
            expected = solver.integrateInnerProduct(z,values,zDomain);

            testCase.verifyEqual(weights.'*values,expected,RelTol=2e-14,AbsTol=2e-14)
        end
    end
end
