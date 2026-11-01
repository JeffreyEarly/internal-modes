classdef IMClassicalRootTransformTests < matlab.unittest.TestCase
    % Independent checks of IMModes.classicalRootTransform (issue #55).

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
        function constantStratificationIsScaledDCTAndDST(testCase)
            N0 = 5e-3; D = 1000;
            background = IMModeBackground.constantStratification(N0=N0,zDomain=[-D 0],surfaceCondition="rigidLid");
            for nPoints = [17 33]
                [transform,grid,assessment] = IMModes.classicalRootTransform(background,nPoints);
                r = nPoints-2;
                dz = D/(nPoints-1);
                % Uniform root points and trapezoid weights.
                testCase.verifyEqual(grid.z,linspace(-D,0,nPoints).',AbsTol=1e-9*D)
                testCase.verifyEqual(grid.weights,dz*[0.5;ones(nPoints-2,1);0.5],RelTol=1e-9)
                testCase.verifyEqual(transform.weights,grid.weights)
                % Full active band, aligned identities and masks.
                testCase.verifyEqual(assessment.retainedInternalModeCount,r)
                testCase.verifyEqual(assessment.limitingCriterion,"maximumRootGridBand")
                testCase.verifyEqual(transform.modeNumber,0:r)
                testCase.verifyEqual(transform.activeModeMask(variable="F"),true(1,r+1))
                testCase.verifyEqual(transform.activeModeMask(variable="G"),[false true(1,r)])
                testCase.verifyLessThan(assessment.final.energyErrorF,1e-10)
                testCase.verifyLessThan(assessment.final.energyErrorG,1e-10)
                % Analytical shapes, equivalent depths and F = h G'.
                g = background.g; z = grid.z; j = 1:r;
                h = N0^2*D^2./(g*(j*pi).^2);
                amplitude = sqrt(2*g/(N0^2*D));
                Gexact = amplitude*sin((z+D)*(j*pi/D));
                Fexact = h.*amplitude.*(j*pi/D).*cos((z+D)*(j*pi/D));
                F = transform.inverseMatrix(variable="F");
                G = transform.inverseMatrix(variable="G");
                testCase.verifyEqual(transform.h(2:end),h,RelTol=1e-10)
                testCase.verifyEqual(F(:,1),ones(nPoints,1))
                testCase.verifyEqual(G(:,1),zeros(nPoints,1))
                % One shared orientation per mode: F and G scale together.
                orientation = sign(sum(G(:,2:end).*Gexact,1));
                testCase.verifyEqual(abs(orientation),ones(1,r))
                testCase.verifyEqual(F(:,2:end),Fexact.*orientation,AbsTol=1e-9*max(abs(Fexact(:))))
                testCase.verifyEqual(G(:,2:end),Gexact.*orientation,AbsTol=1e-9*amplitude)
                testCase.verifyEqual(G([1 end],:),zeros(2,r+1))
            end
        end

        function exponentialStratificationMatchesAnalyticalModes(testCase)
            background = IMModeBackground.exponentialStratification(N0=5.2e-3,b=1300,zDomain=[-4000 0],surfaceCondition="rigidLid");
            [transform,grid,assessment] = IMModes.classicalRootTransform(background,33);
            r = assessment.retainedInternalModeCount;
            solution = background.analyticalSolution();
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) solution.N2(z),zDomain=background.zDomain, ...
                rotationRate=background.rotationRate,latitude=background.latitude,g=background.g);
            exact = solution.internalModes(evp,nModes=r,normalization="geostrophic");
            testCase.verifyEqual(transform.h(2:end),exact.h(1:r),RelTol=1e-8)
            F = transform.inverseMatrix(variable="F"); G = transform.inverseMatrix(variable="G");
            Fexact = exact.F(grid.z,normalization="geostrophic"); Gexact = exact.G(grid.z,normalization="geostrophic");
            orientation = sign(sum(G(:,2:end).*Gexact,1));
            testCase.verifyEqual(F(:,2:end),Fexact.*orientation,AbsTol=1e-6*max(abs(Fexact(:))))
            testCase.verifyEqual(G(:,2:end),Gexact.*orientation,AbsTol=1e-6*max(abs(Gexact(:))))
            testCase.verifyLessThanOrEqual(assessment.final.energyErrorF,assessment.gramTolerance)
            testCase.verifyLessThanOrEqual(assessment.final.energyErrorG,assessment.gramTolerance)
        end

        function rootPointsReproduceV4Grids(testCase)
            fixture = readtable(fullfile(fileparts(mfilename("fullpath")),"Fixtures","IMClassicalRootTransformV4Points.csv"),TextType="string");
            profiles = struct(gmExponential=@(z) (5.2e-3)^2*exp(2*z/1300),pycnocline=@(z) 1e-6+4e-5*exp(-((z+300)/80).^2));
            cases = unique(fixture(:,["profile","depth","nPoints"]),"rows");
            for k = 1:height(cases)
                rows = sortrows(fixture(fixture.profile==cases.profile(k) & fixture.nPoints==cases.nPoints(k),:),"node");
                background = IMModeBackground(N2=profiles.(cases.profile(k)),zDomain=[-cases.depth(k) 0],latitude=33,surfaceCondition="rigidLid");
                [~,grid] = IMModes.classicalRootTransform(background,cases.nPoints(k));
                % v4 roots carry their own unqualified resolution error, so
                % agreement is required to the declared root tolerance.
                testCase.verifyEqual(grid.z,rows.z,cases.profile(k)+" "+cases.nPoints(k),AbsTol=1e-6*cases.depth(k))
            end
        end

        function stronglyVaryingProfileIsIndependentOfResolution(testCase)
            background = IMModeBackground(N2=@(z) 1e-6+4e-5*exp(-((z+300)/80).^2),zDomain=[-2000 0],latitude=33,surfaceCondition="rigidLid");
            [T1,grid1,A1] = IMModes.classicalRootTransform(background,33);
            % The default path raises resolution until an independent refinement agrees.
            testCase.verifyGreaterThan(height(A1.convergence.attempts),1)
            testCase.verifyEqual(A1.convergence.attempts.nEVP(end),A1.convergence.candidateNEVP)
            [T2,grid2,A2] = IMModes.classicalRootTransform(background,33, ...
                solver=IMSolverSpectral(nEVP=2*A1.convergence.candidateNEVP,coordinateKind="wkb"));
            testCase.verifyEqual(height(A2.convergence.attempts),1)
            testCase.verifyEqual(A2.retainedInternalModeCount,A1.retainedInternalModeCount)
            testCase.verifyEqual(grid2.z,grid1.z,AbsTol=1e-8*2000)
            testCase.verifyEqual(T2.h,T1.h,RelTol=1e-8)
            testCase.verifyEqual(grid2.weights,grid1.weights,AbsTol=1e-6*max(grid1.weights))
            testCase.verifyLessThanOrEqual(A1.convergence.maximumShapeError,A1.convergence.tolerance)
            testCase.verifyLessThanOrEqual(A1.convergence.rootDisplacement,A1.convergence.tolerance)
            % A supplied solver is used as given and must pass the check itself.
            testCase.verifyError(@() IMModes.classicalRootTransform(background,33,solver=IMSolverSpectral(nEVP=128,coordinateKind="wkb")), ...
                "IMModes:UnconvergedClassicalModes")
        end

        function automaticSelectionRecordsSearchAndLimit(testCase)
            background = IMModeBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
            [transform,grid,assessment] = IMModes.classicalRootTransform(background,33);
            search = assessment.search;
            r = assessment.retainedInternalModeCount;
            testCase.verifyEqual(assessment.selection,"automatic")
            testCase.verifyEqual(search.internalModeCount.',31:-1:r)
            testCase.verifyEqual(search.accepted.',[false(1,height(search)-1) true])
            testCase.verifyLessThan(r,31)
            testCase.verifyEqual(assessment.limitingCriterion,search.rejectionReason(end-1))
            testCase.verifyEqual(assessment.final,search(end,:))
            testCase.verifyEqual(numel(transform.modeNumber),r+1)
            testCase.verifyEqual(grid.kind,"modeRoots")
            testCase.verifyEqual(grid.design.retainedInternalModeCount,r)
            % The returned weights are the final band's fit, not a full-band fit.
            Phi = transform.inverseMatrix(variable="F");
            gram = Phi.'*(grid.weights.*Phi);
            H = [1000 transform.h(2:end)];
            epsF = norm((gram-diag(H))./sqrt(H.'*H),2);
            testCase.verifyEqual(epsF,assessment.final.energyErrorF,RelTol=1e-8)
        end

        function strictRequestsAreExact(testCase)
            background = IMModeBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
            [transform,~,assessment] = IMModes.classicalRootTransform(background,33,nInternalModes=20);
            testCase.verifyEqual(assessment.selection,"strict")
            testCase.verifyEqual(assessment.retainedInternalModeCount,20)
            testCase.verifyEqual(height(assessment.search),1)
            testCase.verifyEqual(transform.modeNumber,0:20)
            testCase.verifyError(@() IMModes.classicalRootTransform(background,33,nInternalModes=31),"IMModes:ClassicalBandRejected")
            testCase.verifyError(@() IMModes.classicalRootTransform(background,33,nInternalModes=32),"IMModes:InvalidClassicalBand")
        end

        function galerkinRecoveryDiffersFromPrescribedDual(testCase)
            background = IMModeBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
            [transform,grid] = IMModes.classicalRootTransform(background,33);
            Phi = transform.inverseMatrix(variable="F");
            coefficients = linspace(1,-1,size(Phi,2)).';
            values = Phi*coefficients;
            galerkin = transform.transformForward(values,variable="F");
            H = [1000 transform.h(2:end)].';
            prescribed = transform.modeProjectionFunctional(values,variable="F")./H;
            testCase.verifyEqual(galerkin,coefficients,AbsTol=1e-10*norm(coefficients))
            testCase.verifyGreaterThan(norm(prescribed-coefficients),1e3*norm(galerkin-coefficients))
            % A nonmodal input: Galerkin and prescribed-dual coefficients differ.
            x = exp(grid.z/200);
            testCase.verifyGreaterThan(norm(transform.transformForward(x,variable="F")-transform.modeProjectionFunctional(x,variable="F")./H),0)
        end

        function invalidRequestsAreRejected(testCase)
            free = IMModeBackground.constantStratification(N0=5e-3,zDomain=[-1000 0],surfaceCondition="freeSurface");
            testCase.verifyError(@() IMModes.classicalRootTransform(free,17),"IMModes:ClassicalRequiresRigidLid")
            neutral = IMModeBackground(N2=@(z) 1e-5*(z < -200),zDomain=[-1000 0],surfaceCondition="rigidLid");
            testCase.verifyError(@() IMModes.classicalRootTransform(neutral,17),"IMModes:ClassicalRequiresPositiveStratification")
            rigid = IMModeBackground.constantStratification(N0=5e-3,zDomain=[-1000 0],surfaceCondition="rigidLid");
            testCase.verifyError(@() IMModes.classicalRootTransform(rigid,17,solver=IMSolverFiniteDifference(z=linspace(-1000,0,129).')),"IMModes:ClassicalRequiresSpectralSolver")
        end

        function diagonalMetricForwardUsesEquivalentWeightedQR(testCase)
            z = linspace(-1,0,40).';
            Phi = [ones(40,1) cos(pi*(z+1)*(1:9))];
            w = 0.5+rand(RandStream("twister",Seed=3),40,1); w = w/sum(w);
            projection = IMProjection(Phi,diag(w),diag([1 0.5*ones(1,9)]));
            normal = (Phi.'*(w.*Phi))\(Phi.'.*w.');
            testCase.verifyEqual(projection.forwardMatrix,normal,AbsTol=1e-10*norm(normal))
            testCase.verifyEqual(projection.forwardMatrix*Phi,eye(10),AbsTol=1e-12)
        end
    end
end
