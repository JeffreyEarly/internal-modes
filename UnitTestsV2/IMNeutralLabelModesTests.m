classdef IMNeutralLabelModesTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            originalPath = path;
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            testCase.addTeardown(@() path(originalPath));
        end
    end

    methods (Test)
        function analyticalProfilesAndDerivatives(testCase)
            intervals = {[-800 -600],[-200 0],[-1000 -850],[-1000 0]};
            types = ["DD","DN","ND","NN"];
            for i = 1:numel(intervals)
                bounds = intervals{i};
                basis = IMNeutralLabelModesBasis(bounds,zDomain=[-1000 0],nModes=3,constantLabelPolicy="exclude");
                L = diff(bounds);
                n = 1:3;
                if types(i) == "DN" || types(i) == "ND"
                    n = n-0.5;
                end
                m = pi*n/L;
                z = linspace(bounds(1),bounds(2),31).';
                x = z-bounds(1);
                amplitude = sqrt(2000/L);
                if startsWith(types(i),"D")
                    expectedG = amplitude*sin(x*m);
                    expectedGz = amplitude*m.*cos(x*m);
                    expectedGzz = -amplitude*m.^2.*sin(x*m);
                else
                    expectedG = amplitude*cos(x*m);
                    expectedGz = -amplitude*m.*sin(x*m);
                    expectedGzz = -amplitude*m.^2.*cos(x*m);
                end
                testCase.verifyEqual(basis.m_n,m,AbsTol=1e-16)
                testCase.verifyEqual(basis.G(z),expectedG,AbsTol=5e-14)
                testCase.verifyEqual(basis.Gz(z),expectedGz,AbsTol=2e-15)
                testCase.verifyEqual(basis.Gzz(z),expectedGzz,AbsTol=1e-16)
                testCase.verifyEqual(basis.boundaryTypes,[extractBefore(types(i),2),extractAfter(types(i),1)])
                testCase.verifyEqual(basis.APV(z),-basis.f0*expectedGz,AbsTol=1e-18)
            end
        end

        function independentContinuousNormalization(testCase)
            intervals = [-1000 -850;-700 -400;-250 0];
            basis = IMNeutralLabelModesBasis(intervals,zDomain=[-1000 0],nModes=[2 3 4]);
            sampledGGram = zeros(9);
            sampledQGram = zeros(9);
            % Independent adaptive quadrature, one interval at a time.
            for i = 1:9
                for j = 1:9
                    for k = 1:3
                        sampledGGram(i,j) = sampledGGram(i,j) + integral(@(z) profileProduct(basis,z,i,j,"G"),intervals(k,1),intervals(k,2),AbsTol=1e-11)/1000;
                        sampledQGram(i,j) = sampledQGram(i,j) + integral(@(z) profileProduct(basis,z,i,j,"APV"),intervals(k,1),intervals(k,2),AbsTol=1e-20);
                    end
                end
            end
            expectedQGram = diag(basis.f0^2*1000*basis.m_n.^2);
            testCase.verifyEqual(sampledGGram,eye(9),AbsTol=2e-13)
            testCase.verifyEqual(sampledQGram,expectedQGram,AbsTol=2e-20)
            testCase.verifyEqual(basis.gramMatrix(),eye(9))
            testCase.verifyEqual(basis.gramMatrix(variable="APV"),expectedQGram)
        end

        function intervalSupportLabelsAndCollection(testCase)
            basis = IMNeutralLabelModesBasis([-1000 -850;-700 -400;-250 0],zDomain=[-1000 0],nModes=[2 1 3]);
            testCase.verifyEqual(basis.intervalIndex,[1 1 2 3 3 3])
            testCase.verifyEqual(basis.localModeNumber,[1 2 1 1 2 3])
            testCase.verifyEqual(basis.columnLabels,["interval1:mode1","interval1:mode2","interval2:mode1","interval3:mode1","interval3:mode2","interval3:mode3"])
            z = [-1000;-900;-850;-800;-700;-550;-400;-350;-250;-125;0];
            G = basis.G(z);
            testCase.verifyEqual(G([4 8],:),zeros(2,6))
            testCase.verifyEqual(G(1:3,3:6),zeros(3,4))
            testCase.verifyEqual(G(9:11,1:3),zeros(3,3))
            % Interior Gz uses the neutral-side trace; just outside is zero.
            derivative = basis.Gz([-850;-850+1e-6]);
            testCase.verifyGreaterThan(abs(derivative(1,1)),0.01)
            testCase.verifyEqual(derivative(2,1),0)
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            testCase.verifyEqual(collection.metadata.columnKind,"mode")
            testCase.verifyEqual(collection.metadata.columnLabels,basis.columnLabels)
            for derivativeOrder = 0:2
                evaluators = {@() basis.G(z),@() basis.Gz(z),@() basis.Gzz(z)};
                expected = evaluators{derivativeOrder+1}();
                actual = collection.evaluate(z,variable="G",derivativeOrder=derivativeOrder,pages=[2 1 2],columns=[6 2],sampleChunkSize=2,pageChunkSize=1);
                testCase.verifyEqual(actual,repmat(expected(:,[6 2]),1,1,3))
            end
            testCase.verifyEqual(collection.evaluate(z,variable="APV",pages=1),basis.APV(z))
            testCase.verifyError(@() collection.evaluate(z,variable="APV",derivativeOrder=1),"IMBasisCollection:UnsupportedDerivative")
            testCase.verifyError(@() basis.G(-1001),"IMNeutralLabelModesBasis:OutsideDomain")
        end

        function endpointAccelerationPolicy(testCase)
            for acceleration = [0 0.1 -0.1 Inf -Inf]
                basis = IMNeutralLabelModesBasis([-1000 -800;-200 0],zDomain=[-1000 0],nModes=2,g0=acceleration,gd=acceleration);
                G = basis.G([-1000;0]);
                Gz = basis.Gz([-1000;0]);
                if acceleration == 0
                    testCase.verifyEqual(basis.boundaryTypes,["N","D";"D","N"])
                    testCase.verifyEqual(Gz(1,1:2),zeros(1,2),AbsTol=1e-16)
                    testCase.verifyEqual(Gz(2,3:4),zeros(1,2),AbsTol=1e-16)
                    testCase.verifyGreaterThan(min(abs([G(1,1:2),G(2,3:4)])),1)
                else
                    testCase.verifyEqual(basis.boundaryTypes,repmat("D",2,2))
                    testCase.verifyEqual(G,zeros(2,4),AbsTol=2e-14)
                end
            end
            mixed = IMNeutralLabelModesBasis([-1000 -800;-200 0],zDomain=[-1000 0],gd=-Inf,g0=0,nModes=1);
            testCase.verifyEqual(mixed.boundaryTypes,["D","D";"D","N"])
        end

        function constantLabelKernelIsExplicit(testCase)
            testCase.verifyError(@() IMNeutralLabelModesBasis([-1000 0],zDomain=[-1000 0],constantLabelPolicy="reject"),"IMNeutralLabelModesBasis:ConstantLabelKernel")
            basis = IMNeutralLabelModesBasis([-1000 0],zDomain=[-1000 0],nModes=3,constantLabelPolicy="exclude");
            testCase.verifyEqual(basis.m_n,pi*(1:3)/1000,AbsTol=1e-16)
            z = linspace(-1000,0,1001).';
            testCase.verifyEqual(trapz(z,basis.G(z)),zeros(1,3),AbsTol=1e-12)
            testCase.verifyFalse(isprop(basis,"h"))
            testCase.verifyFalse(isprop(basis,"eigenvalues"))
            % A nonzero surface acceleration removes the constant kernel.
            fixed = IMNeutralLabelModesBasis([-1000 0],zDomain=[-1000 0],g0=9.81,nModes=1);
            testCase.verifyEqual(fixed.boundaryTypes,["N","D"])
        end

        function constantLabelIsAutomaticAndDoesNotConsumeCount(testCase)
            basis = IMNeutralLabelModesBasis([-1000 0],zDomain=[-1000 0],nModes=0);
            z = linspace(-1000,0,17).';
            testCase.verifyEqual(basis.constantLabelPolicy,"include")
            testCase.verifyEqual(basis.normalization,"depthRMSG")
            testCase.verifyEqual(basis.nModes,0)
            testCase.verifyEqual(basis.intervalIndex,1)
            testCase.verifyEqual(basis.localModeNumber,0)
            testCase.verifyEqual(basis.columnLabels,"interval1:mode0")
            testCase.verifyEqual(basis.m_n,0)
            testCase.verifyEqual(basis.G(z),ones(size(z)))
            testCase.verifyEqual(basis.Gz(z),zeros(size(z)))
            testCase.verifyEqual(basis.Gzz(z),zeros(size(z)))
            testCase.verifyEqual(basis.APV(z),zeros(size(z)))
            testCase.verifyEqual(basis.gramMatrix(),1)
            testCase.verifyEqual(basis.gramMatrix(variable="APV"),0)
            testCase.verifyError(@() basis.projectionOnGrid(z,ones(size(z))),"IMNeutralLabelModesBasis:UnobservableConstantLabel")
            regular = IMNeutralLabelModesBasis([-500 -400;-200 0],zDomain=[-1000 0],nModes=[0 2]);
            testCase.verifyEqual(regular.intervalIndex,[2 2])
            testCase.verifyEqual(regular.localModeNumber,[1 2])
        end

        function constantLabelPositiveSubsetHasIndependentAPVDual(testCase)
            basis = IMNeutralLabelModesBasis([-1000 0],zDomain=[-1000 0],nModes=3);
            testCase.verifyEqual(basis.localModeNumber,[0 1 2 3])
            testCase.verifyEqual(basis.m_n,pi*(0:3)/1000,AbsTol=1e-16)
            [z,w] = intervalMidpoints([-1000 0],101);
            testCase.verifyError(@() basis.projectionOnGrid(z,w),"IMNeutralLabelModesBasis:UnobservableConstantLabel")
            testCase.verifyError(@() basis.projectionOnGrid(z,w,columns=[2 1]),"IMNeutralLabelModesBasis:UnobservableConstantLabel")
            columns = [4 2 3];
            projection = basis.projectionOnGrid(z,w,columns=columns);
            m = pi*(1:3)/1000;
            q = basis.f0*sqrt(2)*m.*sin((z+1000)*m);
            coefficients = [1+2i;-3i;2];
            testCase.verifyEqual(projection.project(q*coefficients),coefficients(columns-1),AbsTol=2e-13)
            testCase.verifyEqual(projection.targetGramMatrix,eye(3))
            testCase.verifyEqual(projection.majorantGramMatrix,diag(basis.f0^2*1000*basis.m_n(columns).^2))
            testCase.verifyEqual(projection.provenance.unobservableColumns,1)
            testCase.verifyEqual(projection.columnLabels,basis.columnLabels(columns))
            collection = IMBasisCollection({basis});
            testCase.verifyEqual(collection.evaluate(z,variable="G",columns=1),ones(size(z)))
            testCase.verifyError(@() collection.projectionOnGrid(z,w,variable="APV"),"IMNeutralLabelModesBasis:UnobservableConstantLabel")
            selected = collection.projectionOnGrid(z,w,variable="APV",columns=columns);
            testCase.verifyEqual(selected.project(q*coefficients),coefficients(columns-1),AbsTol=2e-13)
        end

        function complexAPVProjectionAcrossIntervals(testCase)
            intervals = [-1000 -850;-700 -400;-250 0];
            basis = IMNeutralLabelModesBasis(intervals,zDomain=[-1000 0],nModes=[2 3 4]);
            [z,w] = intervalMidpoints(intervals,101);
            coefficients = [(1:9).' + 1i*(9:-1:1).',exp(1i*(1:9).')];
            projection = basis.projectionOnGrid(z,w);
            testCase.verifyEqual(projection.project(basis.APV(z)*coefficients),coefficients,AbsTol=1e-12)
            testCase.verifyEqual(projection.projectionKind,"prescribedDual")
            testCase.verifyEqual(projection.targetGramMatrix,eye(9))
            testCase.verifyEqual(projection.majorantGramMatrix,basis.gramMatrix(variable="APV"))
            testCase.verifyEqual(projection.provenance.continuousAPVGramMatrix,basis.gramMatrix(variable="APV"))
            testCase.verifyFalse(projection.supportsGramAssessment)
            testCase.verifyFalse(projection.supportsRoundTripAssessment)
            collection = IMBasisCollection({basis});
            selected = collection.projectionOnGrid(z,w,variable="APV",columns=[8 2]);
            testCase.verifyEqual(selected.project(basis.APV(z)*coefficients),coefficients([8 2],:),AbsTol=1e-12)
            testCase.verifyEqual(selected.columnLabels,basis.columnLabels([8 2]))
            testCase.verifyEqual(collection.projectionOnGrid(z,w).forwardMatrix,projection.forwardMatrix)
            recipe = basis.projectionRecipe();
            testCase.verifyEqual(recipe.targetGramMatrix,projection.targetGramMatrix)
            testCase.verifyFalse(recipe.isAvailable)
            testCase.verifyFalse(recipe.supportsQuadratic)
            testCase.verifyError(@() IMProjection.fromRecipe(recipe,z,w),"IMProjection:UnavailableProjection")
            testCase.verifyEqual(collection.assess(z,w).identity.family,"neutralLabelModes")
        end

        function prescribedDualPreservesQuadratureError(testCase)
            basis = IMNeutralLabelModesBasis([-500 -100],zDomain=[-1000 0],nModes=3);
            z = [-500;-490;-410;-260;-100];
            w = [(z(2)-z(1))/2;(z(3:end)-z(1:end-2))/2;(z(end)-z(end-1))/2];
            projection = basis.projectionOnGrid(z,w);
            m = pi*(1:3)/400;
            q = -basis.f0*sqrt(5)*m.*cos((z+500)*m);
            target = diag(basis.f0^2*1000*m.^2);
            expectedForward = target\(q.'.*w.');
            testCase.verifyEqual(projection.forwardMatrix,expectedForward,RelTol=5e-15)
            coefficients = [1+2i;-0.5i;2];
            actual = projection.project(q*coefficients);
            testCase.verifyEqual(actual,expectedForward*q*coefficients,AbsTol=1e-14)
            testCase.verifyGreaterThan(norm(actual-coefficients),0.1)
        end

        function nonmodalAPVMatchesIndependentQuadrature(testCase)
            % Audit the continuous projector with a source outside the retained span.
            intervals = {[-600 -200],[-200 0],[-1000 -700],[-1000 0],[-200 0],[-1000 -700],[-200 0]};
            accelerations = [0 0;0 0;0 0;0 0;-Inf 0;0 0.1;0 0];
            lowerIsDirichlet = [true true false false true true true];
            halfInteger = [false true true false false false true];
            rotations = [repmat(7.2921e-5,1,6),1e-12];
            % Independent Gauss-Legendre rule from its Jacobi matrix.
            beta = (1:63)./sqrt(4*(1:63).^2-1);
            [V,Lambda] = eig(diag(beta,1)+diag(beta,-1));
            [nodes,order] = sort(diag(Lambda));
            weights = 2*V(1,order).'.^2;
            for i = 1:numel(intervals)
                bounds = intervals{i};
                L = diff(bounds);
                basis = IMNeutralLabelModesBasis(bounds,zDomain=[-1000 0],nModes=7,g0=accelerations(i,1),gd=accelerations(i,2),constantLabelPolicy="exclude",latitude=-45,rotationRate=rotations(i));
                z = bounds(1)+L*(nodes+1)/2;
                f = 2*rotations(i)*sind(-45);
                source = @(x) 1+0.25i+0.7*x/L+0.4*(x/L).^3;
                actual = basis.projectionOnGrid(z,L*weights/2).project(f*source(z-bounds(1)));
                expected = zeros(7,1);
                for j = 1:7
                    m = pi*(j-0.5*halfInteger(i))/L;
                    if lowerIsDirichlet(i)
                        derivative = @(x) sqrt(2000/L)*m*cos(m*x);
                    else
                        derivative = @(x) -sqrt(2000/L)*m*sin(m*x);
                    end
                    % Cancel f analytically; adaptive integration is independent
                    % of both the basis evaluator and the supplied sample rule.
                    expected(j) = -integral(@(x) derivative(x).*source(x),0,L,AbsTol=1e-12,RelTol=1e-12)/(1000*m^2);
                end
                testCase.verifyLessThan(norm(actual-expected)/max(1,norm(expected)),2e-12)
            end
        end

        function constantAPVHasIndependentEndpointIntegral(testCase)
            intervals = {[-800 -600],[-200 0],[-1000 -800]};
            for i = 1:numel(intervals)
                bounds = intervals{i};
                basis = IMNeutralLabelModesBasis(bounds,zDomain=[-1000 0],nModes=3);
                [z,w] = intervalMidpoints(bounds,10001);
                projection = basis.projectionOnGrid(z,w);
                endpointG = basis.G(bounds.');
                expected = -diff(endpointG,1,1).'./(basis.f0*1000*basis.m_n(:).^2);
                actual = projection.project(ones(size(z)));
                testCase.verifyEqual(actual,expected,AbsTol=1e-4,RelTol=1e-7)
                if i == 1
                    testCase.verifyEqual(actual,zeros(3,1),AbsTol=1e-9)
                else
                    testCase.verifyGreaterThan(min(abs(actual)),1e4)
                end
            end
        end

        function hemisphereAndSmallNonzeroRotation(testCase)
            coefficients = [1+1i;2-1i;0.5i];
            [z,w] = intervalMidpoints([-400 -200],101);
            for rotationRate = [7.2921e-5 1e-8]
                north = IMNeutralLabelModesBasis([-400 -200],zDomain=[-1000 0],nModes=3,rotationRate=rotationRate,latitude=33);
                south = IMNeutralLabelModesBasis([-400 -200],zDomain=[-1000 0],nModes=3,rotationRate=rotationRate,latitude=-33);
                testCase.verifyEqual(south.APV(z),-north.APV(z))
                for basis = {north,south}
                    projection = basis{1}.projectionOnGrid(z,w);
                    testCase.verifyEqual(projection.project(basis{1}.APV(z)*coefficients),coefficients,AbsTol=2e-13)
                    testCase.verifyEqual(projection.sampledGramRank,3)
                end
            end
            for options = {struct("rotationRate",0),struct("latitude",0)}
                args = namedargs2cell(options{1});
                basis = IMNeutralLabelModesBasis([-400 -200],args{:},zDomain=[-1000 0],nModes=3);
                testCase.verifyEqual(basis.APV(z),zeros(numel(z),3))
                testCase.verifyEqual(basis.gramMatrix(variable="APV"),zeros(3))
                testCase.verifyError(@() basis.projectionOnGrid(z,w),"IMNeutralLabelModesBasis:ZeroRotationProjection")
            end
        end

        function invalidIntervalsAndProjectionInputs(testCase)
            for intervals = {zeros(0,2),[-100 -200],[-1001 -900],[-200 1],[-500 -200;-400 -100],[-500 -200;-200 -100],[-200 -100;-500 -400]}
                testCase.verifyError(@() IMNeutralLabelModesBasis(intervals{1},zDomain=[-1000 0]),"IMNeutralLabelModesBasis:InvalidIntervals")
            end
            testCase.verifyError(@() IMNeutralLabelModesBasis([-500 -400],zDomain=[0 -1000]),"IMNeutralLabelModesBasis:InvalidDomain")
            testCase.verifyError(@() IMNeutralLabelModesBasis([-500 -400],zDomain=[-1000 0],nModes=[1 2]),"IMNeutralLabelModesBasis:InvalidModeCounts")
            testCase.verifyError(@() IMNeutralLabelModesBasis([-500 -400],zDomain=[-1000 0],g0=NaN),"IMNeutralLabelModesBasis:InvalidEndpointAcceleration")
            basis = IMNeutralLabelModesBasis([-500 -400],zDomain=[-1000 0],nModes=3);
            testCase.verifyError(@() basis.projectionOnGrid([-500;-500],[1;1]),"IMNeutralLabelModesBasis:InvalidGrid")
            testCase.verifyError(@() basis.projectionOnGrid([-500;-400],1),"IMNeutralLabelModesBasis:InvalidWeights")
            testCase.verifyError(@() basis.projectionOnGrid([-500;-400],[0;0]),"IMNeutralLabelModesBasis:InvalidWeights")
            testCase.verifyError(@() basis.projectionOnGrid([-500;-400],[1;1],columns=[1 1]),"IMProjection:InvalidColumns")
            testCase.verifyError(@() basis.projectionOnGrid([-500;-400],[1;1],columns=4),"IMProjection:InvalidColumns")
            testCase.verifyError(@() basis.projectionOnGrid([-1001;0],[1;1]),"IMNeutralLabelModesBasis:OutsideDomain")
        end
    end
end

function product = profileProduct(basis,z,i,j,variable)
values = basis.(variable)(z(:));
product = reshape(values(:,i).*values(:,j),size(z));
end

function [z,w] = intervalMidpoints(intervals,count)
z = zeros(count*size(intervals,1),1);
w = z;
for i = 1:size(intervals,1)
    rows = (i-1)*count+(1:count);
    dz = diff(intervals(i,:))/count;
    z(rows) = intervals(i,1) + ((1:count).'-0.5)*dz;
    w(rows) = dz;
end
end
