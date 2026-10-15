classdef IMModeBackgroundTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function arbitraryBackgroundOwnsThePhysicalColumn(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            background = IMModeBackground(N2=N2,zDomain=[-1000 0]);
            testCase.verifyEqual(background.N2,N2);
            testCase.verifyEqual(background.zDomain,[-1000 0]);
            testCase.verifyEqual(background.D,1000);
            testCase.verifyEqual(background.rotationRate,7.2921e-5);
            testCase.verifyEqual(background.latitude,33);
            testCase.verifyEqual(background.f0,2*7.2921e-5*sind(33));
            testCase.verifyEqual(background.g,9.81);
            testCase.verifyEqual(background.surfaceCondition,"freeSurface");
            testCase.verifyEqual(background.profileKind,"arbitrary");
            testCase.verifyEqual(background.profileParameters,struct());
            testCase.verifyError(@() background.analyticalSolution(),"IMModeBackground:UnsupportedAnalyticalProfile");
        end

        function sharedPhysicsPreservesHemisphereAndDomainOrigin(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[100 500],rotationRate=1e-4,latitude=-30,g=10,surfaceCondition="rigidLid");
            testCase.verifyEqual(background.D,400);
            testCase.verifyEqual(background.f0,-1e-4,AbsTol=1e-19);
            testCase.verifyEqual(background.g,10);
            testCase.verifyEqual(background.surfaceCondition,"rigidLid");
        end

        function constantSourceCreatesMatchingAnalyticalPhysics(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 100],rotationRate=1e-4,latitude=-20,g=10,surfaceCondition="rigidLid");
            z = reshape(linspace(-1000,100,6),[2 3]);
            testCase.verifyEqual(background.N2(z),1e-4*ones(size(z)));
            testCase.verifyEqual(background.profileKind,"constant");
            testCase.verifyEqual(background.profileParameters,struct("N0",.01));
            solution = background.analyticalSolution();
            testCase.verifyClass(solution,"IMConstantStratificationSolution");
            testCase.verifyEqual(solution.N2(z),background.N2(z));
            testCase.verifyEqual(solution.zDomain,background.zDomain);
            testCase.verifyEqual(solution.rotationRate,background.rotationRate);
            testCase.verifyEqual(solution.latitude,background.latitude);
            testCase.verifyEqual(solution.g,background.g);
        end

        function exponentialSourceAnchorsN0AtTheSurface(testCase)
            background = IMModeBackground.exponentialStratification(N0=.01,b=500,zDomain=[-900 100]);
            z = [-900 -400 100];
            testCase.verifyEqual(background.N2(z),1e-4*exp(2*(z-100)/500),RelTol=1e-14);
            testCase.verifyEqual(background.N2(100),1e-4);
            testCase.verifyEqual(background.profileKind,"exponential");
            testCase.verifyEqual(background.profileParameters,struct("N0",.01,"b",500));
            testCase.verifyError(@() background.analyticalSolution(),"IMModeBackground:UnsupportedAnalyticalDomain");
        end

        function exponentialAnalyticalSourceHasOneProfile(testCase)
            background = IMModeBackground.exponentialStratification(N0=.02,b=700,zDomain=[-1000 0],rotationRate=0,latitude=0,g=10);
            solution = background.analyticalSolution();
            z = linspace(-1000,0,11).';
            testCase.verifyClass(solution,"IMExponentialStratificationSolution");
            testCase.verifyEqual(solution.N2(z),background.N2(z),RelTol=1e-14);
            testCase.verifyEqual(solution.N0,background.profileParameters.N0);
            testCase.verifyEqual(solution.b,background.profileParameters.b);
            testCase.verifyEqual(solution.zDomain,background.zDomain);
            testCase.verifyEqual(solution.rotationRate,background.rotationRate);
            testCase.verifyEqual(solution.latitude,background.latitude);
            testCase.verifyEqual(solution.g,background.g);
        end

        function domainAndSurfaceVocabularyAreValidated(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            testCase.verifyError(@() IMModeBackground(N2=N2,zDomain=[0 -1000]),"IMModeBackground:InvalidDomain");
            testCase.verifyError(@() IMModeBackground(N2=N2,zDomain=[0 0]),"IMModeBackground:InvalidDomain");
            testCase.verifyError(@() IMModeBackground(N2=N2,zDomain=[-1 Inf]),"MATLAB:validators:mustBeFinite");
            testCase.verifyError(@() IMModeBackground(N2=N2,zDomain=[-1 0],surfaceCondition="free"),"IMModeBackground:InvalidSurfaceCondition");
        end

        function physicalAndSourceParametersAreValidated(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            testCase.verifyError(@() IMModeBackground(N2=N2,zDomain=[-1 0],rotationRate=-1),"MATLAB:validators:mustBeNonnegative");
            testCase.verifyError(@() IMModeBackground(N2=N2,zDomain=[-1 0],latitude=91),"MATLAB:validators:mustBeLessThanOrEqual");
            testCase.verifyError(@() IMModeBackground(N2=N2,zDomain=[-1 0],g=0),"MATLAB:validators:mustBePositive");
            testCase.verifyError(@() IMModeBackground.constantStratification(N0=0),"MATLAB:validators:mustBePositive");
            testCase.verifyError(@() IMModeBackground.exponentialStratification(b=-1),"MATLAB:validators:mustBePositive");
        end

        function energyBoundaryIncludesEveryFiniteEndpoint(testCase)
            both = IMGeneralizedEnergyBoundary(g0=-.02,gd=0);
            testCase.verifyEqual(both.g0,-.02);
            testCase.verifyEqual(both.gd,0);
            testCase.verifyEqual(both.activeEndpoints,["surface", "bottom"]);
            surface = IMGeneralizedEnergyBoundary(g0=0,gd=Inf);
            bottom = IMGeneralizedEnergyBoundary(g0=Inf,gd=.02);
            neither = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf);
            testCase.verifyEqual(surface.activeEndpoints,"surface");
            testCase.verifyEqual(bottom.activeEndpoints,"bottom");
            testCase.verifyEmpty(neither.activeEndpoints);
        end

        function energyBoundaryRejectsNonphysicalInfinityAndNaN(testCase)
            testCase.verifyError(@() IMGeneralizedEnergyBoundary(g0=NaN,gd=0),"IMGeneralizedEnergyBoundary:InvalidCoefficient");
            testCase.verifyError(@() IMGeneralizedEnergyBoundary(g0=0,gd=-Inf),"IMGeneralizedEnergyBoundary:InvalidCoefficient");
        end

        function analyticalSolverDoesNotOwnAnotherProfile(testCase)
            solver = IMSolverAnalytical();
            testCase.verifyClass(solver,"IMSolverAnalytical");
            testCase.verifyEmpty(properties(solver));
        end
    end
end
