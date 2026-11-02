classdef IMZeroAPVCommonBasisTests < matlab.unittest.TestCase
    % Independently check physical forms and common endpoint coordinates.

    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function endpointInventoriesMatchIndependentPhysicalEnergy(testCase)
            coefficients = {[-3 .4],[-3 Inf],[Inf .4],[Inf Inf]};
            for surface = ["freeSurface","rigidLid"]
                background = physicalBackground(surface);
                for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48,coordinateKind="z")}
                    for j = 1:numel(coefficients)
                        boundary = IMGeneralizedEnergyBoundary(g0=coefficients{j}(1),gd=coefficients{j}(2));
                        basis = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=solver{1},energyBoundary=boundary);
                        [H,B0,Bd,response] = independentForms(basis);
                        expectedEndpoints = boundary.activeEndpoints;
                        testCase.verifyEqual(basis.canonicalEndpointLabels,expectedEndpoints);
                        testCase.verifySize(basis.sourceCoordinateMap,[numel(expectedEndpoints),basis.nColumns]);
                        testCase.verifyEqual(basis.energyMatrix,H,RelTol=3e-11,AbsTol=5e-13);
                        testCase.verifyEqual(basis.surfaceBuoyancyMatrix,B0,RelTol=3e-11,AbsTol=5e-13);
                        testCase.verifyEqual(basis.bottomBuoyancyMatrix,Bd,RelTol=3e-11,AbsTol=5e-13);
                        selectedRows = endpointRows(expectedEndpoints);
                        testCase.verifyEqual(response(selectedRows,:),eye(basis.nColumns),AbsTol=3e-12);
                        inactiveRows = setdiff(1:2,selectedRows);
                        testCase.verifyEqual(response(inactiveRows,:),zeros(numel(inactiveRows),basis.nColumns),AbsTol=3e-12);
                        testCase.verifyEqual(basis.endpointResponseMetric,response'*response,AbsTol=3e-12);
                        expectedGeneralized = H;
                        if isfinite(boundary.g0), expectedGeneralized = expectedGeneralized+boundary.g0*B0; end
                        if isfinite(boundary.gd), expectedGeneralized = expectedGeneralized+boundary.gd*Bd; end
                        testCase.verifyEqual(basis.generalizedEnergyMatrix(),expectedGeneralized,RelTol=3e-11,AbsTol=5e-13);
                        testCase.verifyTrue(all(isfinite(basis.generalizedEnergyMatrix()),"all"));
                        if ~ismember("surface",expectedEndpoints)
                            testCase.verifyEqual(basis.surfaceBuoyancyMatrix,zeros(basis.nColumns),AbsTol=0);
                            testCase.verifyEqual(basis.generalizedEnergyMatrix(g0=-1e6),basis.generalizedEnergyMatrix(g0=0),AbsTol=0);
                        end
                        if ~ismember("bottom",expectedEndpoints)
                            testCase.verifyEqual(basis.bottomBuoyancyMatrix,zeros(basis.nColumns),AbsTol=0);
                            testCase.verifyEqual(basis.generalizedEnergyMatrix(gd=-1e6),basis.generalizedEnergyMatrix(gd=0),AbsTol=0);
                        end
                        z = [-1;-.7;-.2;0];
                        testCase.verifyEqual(basis.evaluate(z,"Fz"),-background.N2(z)/background.g.*basis.G(z),AbsTol=3e-12);
                        testCase.verifyEqual(basis.Gz(z),-background.g*basis.k^2/background.f0^2*basis.F(z),AbsTol=3e-12);
                        testCase.verifyEqual(basis.Gzz(z),basis.k^2/background.f0^2*background.N2(z).*basis.G(z),AbsTol=3e-12);
                        testCase.verifyError(@() basis.h,"IMModeBasis:UnavailableEquivalentDepth");
                        testCase.verifyError(@() basis.eigenvalues,"IMModeBasis:UnavailableEigenvalues");
                        testCase.verifyFalse(any(ismember(["h","eigenvalues"],basis.capabilities.spectrumNames)));
                        testCase.verifyError(@() basis.withNormalization("rotationRightMetric"),"IMModeBasis:UnavailableNormalization");
                        if basis.nColumns <= 1
                            rotated = basis.rotateBoundaryDepth();
                            verifyCoordinateChange(testCase,basis,rotated,rotated.rotationMetadata.coordinateChange);
                            testCase.verifyEqual(rotated.nColumns,basis.nColumns);
                        end
                    end
                end
            end
        end

        function boundaryDepthRotationUsesTheSignedManuscriptPencil(testCase)
            background = physicalBackground("freeSurface");
            energy = IMGeneralizedEnergyBoundary(g0=-3,gd=.4);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48,coordinateKind="z")}
                original = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=solver{1},energyBoundary=energy);
                changed = original.rotateBoundaryDepth();
                T = changed.rotationMetadata.coordinateChange;
                lambda = reshape(changed.rotationMetadata.eigenvalues,[],1);
                verifyCoordinateChange(testCase,original,changed,T);
                testCase.verifyTrue(any(lambda < 0));
                testCase.verifyEqual(changed.endpointResponseMetric,eye(2),AbsTol=3e-12);
                testCase.verifyEqual(changed.generalizedEnergyMatrix(),diag(lambda),AbsTol=3e-12);
                testCase.verifyEqual(reshape(changed.rotationMetadata.h0,[],1),2*original.k^2*lambda,AbsTol=3e-13);
                testCase.verifyEqual(reshape(changed.rotationMetadata.rightMetricSignatures,[],1),ones(2,1));
                testCase.verifyEqual(reshape(changed.rotationMetadata.leftMetricSignatures,[],1),sign(lambda));
                testCase.verifyEqual(changed.normalization,"rotationRightMetric");
                testCase.verifyEqual(changed.canonicalEndpointLabels,["surface","bottom"]);
                testCase.verifyNotEqual(changed.columnLabels,changed.canonicalEndpointLabels);
                testCase.verifyError(@() changed.h,"IMModeBasis:UnavailableEquivalentDepth");
                testCase.verifyError(@() changed.eigenvalues,"IMModeBasis:UnavailableEigenvalues");
                [physicalH,B0,Bd] = independentForms(changed);
                testCase.verifyEqual(changed.energyMatrix,physicalH,RelTol=3e-11,AbsTol=3e-12);
                testCase.verifyEqual(changed.generalizedEnergyMatrix(),physicalH+energy.g0*B0+energy.gd*Bd,RelTol=3e-11,AbsTol=3e-12);
            end
        end

        function surfaceAndCustomPencilsUseCurrentCoordinates(testCase)
            background = physicalBackground("freeSurface");
            energy = IMGeneralizedEnergyBoundary(g0=-3,gd=.4);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=48,coordinateKind="z")}
                canonical = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=solver{1},energyBoundary=energy);
                original = canonical.withNormalization("maxAbsF").rotateBoundaryDepth();
                changed = original.rotateSurfaceBuoyancy();
                T = changed.rotationMetadata.coordinateChange;
                lambda = reshape(changed.rotationMetadata.eigenvalues,[],1);
                signatures = reshape(changed.rotationMetadata.rightMetricSignatures,[],1);
                L = background.g*original.surfaceBuoyancyMatrix;
                R = original.generalizedEnergyMatrix();
                testCase.verifyEqual(L*T,R*T*diag(lambda),RelTol=3e-11,AbsTol=3e-12);
                testCase.verifyEqual(T'*R*T,diag(signatures),AbsTol=3e-12);
                testCase.verifyEqual(reshape(changed.rotationMetadata.leftMetricSignatures,[],1),sign(lambda).*signatures);
                testCase.verifyEmpty(changed.rotationMetadata.h0);
                verifyCoordinateChange(testCase,original,changed,T);
                left = [4 .7;.7 -1]; right = [2 .2;.2 1];
                custom = changed.rotateWithPencil(name="currentControl",leftMatrix=left,rightMatrix=right);
                U = custom.rotationMetadata.coordinateChange;
                eigenvalues = reshape(custom.rotationMetadata.eigenvalues,[],1);
                testCase.verifyEqual(left*U,right*U*diag(eigenvalues),RelTol=3e-11,AbsTol=3e-12);
                testCase.verifyEqual(U'*right*U,eye(2),AbsTol=3e-12);
                verifyCoordinateChange(testCase,changed,custom,U);
                named = changed.rotateWithPencil(name="boundaryDepth",leftMatrix=left,rightMatrix=right);
                testCase.verifyEmpty(named.rotationMetadata.h0);
            end
        end

        function selectionAndNormalizationKeepOnlyTheRetainedSubspace(testCase)
            background = physicalBackground("freeSurface");
            energy = IMGeneralizedEnergyBoundary(g0=-3,gd=.4);
            original = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=IMSolverAnalytical(),energyBoundary=energy);
            rotated = original.rotateBoundaryDepth();
            profiled = rotated.withNormalization("maxAbsG");
            restored = profiled.withNormalization("rotationRightMetric");
            testCase.verifyEqual(restored.F([-1;-.4;0]),rotated.F([-1;-.4;0]),RelTol=3e-12,AbsTol=3e-12);
            selected = profiled.selectColumns(2);
            changed = selected.rotateWithPencil(name="oneDirection",leftMatrix=3,rightMatrix=2);
            verifyCoordinateChange(testCase,selected,changed,changed.rotationMetadata.coordinateChange);
            testCase.verifyEqual(changed.nColumns,1);
            testCase.verifySize(changed.sourceCoordinateMap,[2 1]);
            testCase.verifyEqual(rank([selected.sourceCoordinateMap,changed.sourceCoordinateMap]),1);
            testCase.verifyEqual(changed.canonicalEndpointLabels,["surface","bottom"]);
            normalized = changed.withNormalization("maxAbsF");
            ratio = diag(changed.scientificSpecification.resolved.normalization.scales./normalized.scientificSpecification.resolved.normalization.scales);
            verifyCoordinateChange(testCase,changed,normalized,ratio);
            empty = normalized.selectColumns([]);
            emptyRotation = empty.rotateBoundaryDepth();
            testCase.verifyEqual(emptyRotation.nColumns,0);
            testCase.verifySize(emptyRotation.sourceCoordinateMap,[2 0]);
            testCase.verifyEqual(emptyRotation.energyMatrix,zeros(0),AbsTol=0);
            testCase.verifyEqual(emptyRotation.canonicalEndpointLabels,["surface","bottom"]);
            testCase.verifySize(emptyRotation.G([-1;0]),[2 0]);
            testCase.verifyError(@() emptyRotation.withNormalization("endpointResponse"),"IMModeBasis:UnavailableNormalization");
        end

        function collectionsDistinguishEndpointIdentityFromRotatedColumns(testCase)
            background = physicalBackground("freeSurface");
            original = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=IMSolverAnalytical(),energyBoundary=IMGeneralizedEnergyBoundary(g0=-3,gd=.4));
            rotated = original.rotateBoundaryDepth().withNormalization("maxAbsF");
            for basis = {rotated,rotated.selectColumns(2),rotated.selectColumns([])}
                current = basis{1};
                collection = IMBasisCollection({current},basisIndex=[1 1 1]);
                metadata = collection.metadata;
                testCase.verifyEqual(metadata.columnKind,"endpoint");
                testCase.verifyEqual(metadata.endpointLabels,current.canonicalEndpointLabels);
                testCase.verifyEqual(metadata.columnLabels,current.columnLabels);
                testCase.verifyEqual(metadata.sourceCoordinateMap,current.sourceCoordinateMap,AbsTol=0);
                testCase.verifyEqual(metadata.rotationMetadata,current.rotationMetadata);
                testCase.verifyEqual(collection.kappa,[1.2 1.2 1.2]);
                z = [-1;-.7;-.1;0];
                for variable = current.capabilities.variables
                    j = find(current.capabilities.variables == variable,1);
                    for order = current.capabilities.derivativeOrders{j}
                        expected = current.evaluate(z,variable,order);
                        actual = collection.evaluate(z,variable=variable,derivativeOrder=order,pages=[3 1 3 2],sampleChunkSize=2,pageChunkSize=2);
                        testCase.verifyEqual(actual,repmat(expected,1,1,4),AbsTol=3e-12);
                    end
                end
            end
        end

        function physicalDualAndSharedRecoveryCovaryForComplexOffSpanData(testCase)
            background = physicalBackground("freeSurface");
            energy = IMGeneralizedEnergyBoundary(g0=-3,gd=.4);
            original = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=IMSolverAnalytical(),energyBoundary=energy).withNormalization("maxAbsF");
            changed = original.rotateSurfaceBuoyancy();
            T = changed.rotationMetadata.coordinateChange;
            [z,w] = legpts(48,background.zDomain); w = w(:);
            [A,metric] = physicalObservation(original,z,w,energy);
            [newA,newMetric] = physicalObservation(changed,z,w,energy);
            Hplus = A'*(abs(metric).*A);
            U = diag(exp(1i*linspace(.1,1.8,size(A,1))));
            A = U*A; newA = U*newA;
            W = diag(metric);
            H = original.generalizedEnergyMatrix();
            newH = changed.generalizedEnergyMatrix();
            B = H\(A'*W); newB = newH\(newA'*W);
            C = chol(Hplus,"lower")'\eye(original.nColumns);
            D = diag(sqrt(abs(metric)));
            F = diag(linspace(.7,1.6,numel(metric)))*D;
            complement = null(A');
            source = A*[1+.3i;-.2+1i]+(.4-.7i)*complement(:,1)+(.3+.2i)*complement(:,2);
            testCase.verifyEqual(newMetric,metric,AbsTol=0);
            testCase.verifyEqual(newA,A*T,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(A'*W*A,H,RelTol=3e-11,AbsTol=3e-12);
            testCase.verifyEqual(newH,T'*H*T,RelTol=3e-11,AbsTol=3e-12);
            testCase.verifyTrue(any(eig(H) < 0));
            testCase.verifyEqual(newB,T\B,RelTol=3e-11,AbsTol=3e-12);
            testCase.verifyGreaterThan(norm(source-A*(A\source)),.5);
            for policy = ["sampledLeastSquares","correctedPhysicalDual","extendedSpan"]
                oldOptions = struct(policy=policy,columnLabels=original.columnLabels,physicalDual=B,coefficientScale=C,observationFactor=D,fitFactor=F);
                newOptions = struct(policy=policy,columnLabels=changed.columnLabels,physicalDual=newB,coefficientScale=T\C,observationFactor=D,fitFactor=F);
                if policy == "extendedSpan"
                    oldOptions.extendedSynthesis = [A,complement(:,1)];
                    oldOptions.extendedColumnLabels = [original.columnLabels,"latent:1"];
                    oldOptions.extendedActiveColumnMask = true(1,3);
                    oldOptions.extendedCoefficientScale = blkdiag(C,1);
                    oldOptions.retainedColumnIndices = 1:2;
                    newOptions.extendedSynthesis = [newA,complement(:,1)];
                    newOptions.extendedColumnLabels = [changed.columnLabels,"latent:1"];
                    newOptions.extendedActiveColumnMask = true(1,3);
                    newOptions.extendedCoefficientScale = blkdiag(T\C,1);
                    newOptions.retainedColumnIndices = 1:2;
                    testCase.verifyGreaterThan(norm(source-oldOptions.extendedSynthesis*(oldOptions.extendedSynthesis\source)),.3);
                end
                oldPairs = namedargs2cell(oldOptions); newPairs = namedargs2cell(newOptions);
                oldRecovery = IMDiscreteTransform.fromObservationMatrices(A,oldPairs{:});
                newRecovery = IMDiscreteTransform.fromObservationMatrices(newA,newPairs{:});
                testCase.verifyEqual(newRecovery.forwardMatrix,T\oldRecovery.forwardMatrix,RelTol=3e-11,AbsTol=3e-12);
                testCase.verifyEqual(newRecovery.physicalDualMatrix*source,T\(oldRecovery.physicalDualMatrix*source),RelTol=3e-11,AbsTol=3e-12);
                testCase.verifyEqual(newRecovery.transformBack(newRecovery.transformForward(source)),oldRecovery.transformBack(oldRecovery.transformForward(source)),RelTol=3e-11,AbsTol=3e-12);
                testCase.verifyEqual(newRecovery.recoverySpecification.fitFactor,oldRecovery.recoverySpecification.fitFactor,AbsTol=0);
                testCase.verifyEqual(newRecovery.recoverySpecification.observationFactor,oldRecovery.recoverySpecification.observationFactor,AbsTol=0);
                testCase.verifyEqual(newRecovery.diagnostics.observationNoiseGain,oldRecovery.diagnostics.observationNoiseGain,RelTol=3e-11,AbsTol=3e-12);
                if policy == "extendedSpan"
                    testCase.verifyEqual(newRecovery.forwardMatrix*complement(:,1),zeros(2,1),AbsTol=3e-12);
                end
            end
        end

        function transformationsReuseTheComputedEndpointResponses(testCase)
            background = physicalBackground("freeSurface");
            counter = IMGeostrophicZeroAPVTestCounter();
            solver = IMInstrumentedGeostrophicZeroAPVSolver(nEVP=48,solveCounter=counter);
            original = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=solver,energyBoundary=IMGeneralizedEnergyBoundary(g0=-3,gd=.4));
            initialCalls = counter.callCount;
            testCase.verifyEqual(initialCalls,1);
            changed = original.rotateBoundaryDepth().withNormalization("maxAbsG").selectColumns(2).rotateWithPencil(name="retained",leftMatrix=2,rightMatrix=3);
            changed.F([-.9;-.1]); changed.generalizedEnergyMatrix();
            collection = IMBasisCollection({changed},basisIndex=[1 1]);
            collection.evaluate([-.9;-.1],variable="G",pages=[2 1 2]);
            testCase.verifyEqual(counter.callCount,initialCalls);
            testCase.verifyEqual(counter.rightHandSideCounts,2);
        end

        function alreadyRotatedSourcesRetainTheirCanonicalCoordinateMap(testCase)
            background = physicalBackground("freeSurface");
            energy = IMGeneralizedEnergyBoundary(g0=-3,gd=.4);
            source = background.analyticalSolution().geostrophicZeroAPVModesAtWavenumber(1.2).rotateBoundaryDepth(g0=energy.g0,gd=energy.gd);
            info = table(["provided:1";"provided:2"],repmat("endpointResponse",2,1),(1:2).',VariableNames=["columnLabel","branch","branchIndex"]);
            basis = IMModeBasis(source=source,background=background,family="geostrophicZeroAPV",modeInfo=info,normalization="depthRMSF",horizontalKind="fixedWavenumber",k=1.2,frequency=[0 0],hasFrequency=true,scientificInput=struct(energyBoundary=energy));
            scales = basis.scientificSpecification.resolved.normalization.scales;
            z = [-1;-.6;-.1;0];
            testCase.verifyEqual(basis.sourceCoordinateMap,source.rotationMatrix./scales,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(basis.F(z),source.F(z)./scales,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(basis.G(z),source.G(z)./scales,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(basis.energyMatrix,source.energyMatrix./(scales(:)*scales),RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(basis.generalizedEnergyMatrix(),source.generalizedEnergyMatrix(g0=energy.g0,gd=energy.gd)./(scales(:)*scales),RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(basis.canonicalEndpointLabels,["surface","bottom"]);
            testCase.verifyEqual(endpointResponse(basis),basis.sourceCoordinateMap,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyError(@() basis.withNormalization("endpointResponse"),"IMModeBasis:UnavailableNormalization");
            testCase.verifyError(@() IMModeBasis(source=source,background=background,family="geostrophicZeroAPV",modeInfo=info,normalization="depthRMSF",horizontalKind="fixedWavenumber",k=1.3),"IMModeBasis:IncompatibleWavenumber");
            testCase.verifyError(@() IMModeBasis(source=source,background=background,family="geostrophicZeroAPV",modeInfo=info,normalization="depthRMSF",horizontalKind="independent"),"IMModeBasis:IncompatibleWavenumber");
            testCase.verifyError(@() IMModeBasis(background=background,family="geostrophicZeroAPV",normalization="depthRMSF",horizontalKind="fixedWavenumber",k=0),"IMModeBasis:InvalidHorizontalIdentity");
        end
    end
end

function background = physicalBackground(surface)
background = IMModeBackground.constantStratification(N0=1,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,surfaceCondition=surface);
end

function rows = endpointRows(labels)
rows = zeros(1,numel(labels));
rows(labels == "surface") = 1;
rows(labels == "bottom") = 2;
end

function response = endpointResponse(basis)
F = basis.F(basis.zDomain([2 1]).');
G = basis.G(basis.zDomain([2 1]).');
response = G;
if basis.background.surfaceCondition == "freeSurface", response(1,:) = response(1,:)-F(1,:); end
end

function [H,B0,Bd,response] = independentForms(basis)
[z,w] = legpts(64,basis.zDomain); w = w(:);
F = basis.F(z); G = basis.G(z);
f0 = basis.background.f0; g = basis.background.g; k = basis.k;
H = (k^2*(F'*(w.*F))+(f0^2/g^2)*(G'*(w.*basis.background.N2(z).*G)))/(2*k^4);
if basis.background.surfaceCondition == "freeSurface"
    Fs = basis.F(basis.zDomain(2));
    H = H+(f0^2/(2*g*k^4))*(Fs'*Fs);
end
response = endpointResponse(basis);
B0 = (f0^2/(2*g^2*k^4))*(response(1,:)'*response(1,:));
Bd = (f0^2/(2*g^2*k^4))*(response(2,:)'*response(2,:));
end

function verifyCoordinateChange(testCase,original,changed,T)
testCase.verifyEqual(changed.sourceCoordinateMap,original.sourceCoordinateMap*T,RelTol=3e-12,AbsTol=3e-12);
testCase.verifyEqual(changed.canonicalEndpointLabels,original.canonicalEndpointLabels);
testCase.verifyEqual(changed.endpointResponseMetric,T'*original.endpointResponseMetric*T,RelTol=3e-11,AbsTol=3e-12);
testCase.verifyEqual(changed.energyMatrix,T'*original.energyMatrix*T,RelTol=3e-11,AbsTol=3e-12);
testCase.verifyEqual(changed.surfaceBuoyancyMatrix,T'*original.surfaceBuoyancyMatrix*T,RelTol=3e-11,AbsTol=3e-12);
testCase.verifyEqual(changed.bottomBuoyancyMatrix,T'*original.bottomBuoyancyMatrix*T,RelTol=3e-11,AbsTol=3e-12);
testCase.verifyEqual(changed.generalizedEnergyMatrix(),T'*original.generalizedEnergyMatrix()*T,RelTol=3e-11,AbsTol=3e-12);
testCase.verifyEqual(endpointResponse(changed),endpointResponse(original)*T,RelTol=3e-12,AbsTol=3e-12);
z = [-1;-.6;-.1;0];
for variable = original.capabilities.variables
    j = find(original.capabilities.variables == variable,1);
    for order = original.capabilities.derivativeOrders{j}
        testCase.verifyEqual(changed.evaluate(z,variable,order),original.evaluate(z,variable,order)*T,RelTol=3e-12,AbsTol=3e-12);
    end
end
end

function [A,metric] = physicalObservation(basis,z,w,energy)
F = basis.F(z); G = basis.G(z);
Fs = basis.F(basis.zDomain(2)); response = endpointResponse(basis);
f0 = basis.background.f0; g = basis.background.g; k = basis.k;
A = [F;G;Fs;response];
metric = [w/(2*k^2);w.*basis.background.N2(z)*f0^2/(2*g^2*k^4);f0^2/(2*g*k^4);energy.g0*f0^2/(2*g^2*k^4);energy.gd*f0^2/(2*g^2*k^4)];
end
