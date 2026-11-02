classdef IMZeroAPVCoordinateEdgeTests < matlab.unittest.TestCase
    % Check endpoint coordinate provenance and structural admissibility edges.

    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function aPurePermutationDoesNotInventInitialColumnProvenance(testCase)
            background = edgeBackground("freeSurface");
            original = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=IMSolverAnalytical(),energyBoundary=IMGeneralizedEnergyBoundary(g0=0,gd=0));
            changed = original.rotateWithPencil(name="permutation",leftMatrix=diag([1 3]),rightMatrix=eye(2));
            expected = [0 1;1 0];
            testCase.verifyEqual(changed.sourceCoordinateMap,expected,AbsTol=3e-14);
            testCase.verifyEqual(changed.scientificSpecification.resolved.sourceColumns,[2 1]);
            testCase.verifyEqual(changed.constructionRecord.selection.activeSourceColumns,[2 1]);
            testCase.verifyEmpty(changed.constructionRecord.selection.initialColumnIndices);
            testCase.verifyEqual(changed.constructionRecord.selection.sourceDirections,expected,AbsTol=3e-14);
            testCase.verifyEqual(changed.F([-1;-.3;0]),original.F([-1;-.3;0])*expected,AbsTol=3e-13);
            normalized = changed.withNormalization("depthRMSF").selectColumns(1);
            scales = normalized.scientificSpecification.resolved.normalization.scales;
            testCase.verifyEqual(normalized.sourceCoordinateMap,expected(:,1)/scales,AbsTol=3e-13);
            testCase.verifyEqual(normalized.constructionRecord.selection.sourceDirections,expected(:,1),AbsTol=3e-14);
            testCase.verifyEqual(normalized.scientificSpecification.resolved.sourceColumns,2);
            testCase.verifyEmpty(normalized.constructionRecord.selection.initialColumnIndices);
        end

        function selectedDirectSourcesNarrowInheritedRotationEvidence(testCase)
            background = edgeBackground("freeSurface");
            canonical = background.analyticalSolution().geostrophicZeroAPVModesAtWavenumber(1.2);
            sources = {canonical.rotateBoundaryDepth(g0=-3,gd=.4),canonical.rotateWithPencil(name="boundaryNormalized",leftMatrix=[3 .4;.4 -1],rightMatrix=eye(2))};
            for selectedSource = sources
                source = selectedSource{1};
                legacy = IMBasisCollection({source});
                testCase.verifyNotEqual(legacy.metadata.columnLabels,source.endpoints);
                testCase.verifyEqual(legacy.metadata.endpointLabels,source.endpoints);
                testCase.verifyEqual(legacy.metadata.sourceCoordinateMap,source.rotationMatrix,AbsTol=0);
                for columns = {2,[2 1]}
                    chosen = columns{1};
                    basis = wrapSource(source,background,chosen);
                    scales = basis.scientificSpecification.resolved.normalization.scales;
                    expectedMap = source.rotationMatrix(:,chosen)./scales;
                    testCase.verifyEqual(basis.sourceCoordinateMap,expectedMap,AbsTol=3e-12);
                    testCase.verifyTrue(basis.rotationMetadata.isRotated);
                    testCase.verifyEqual(basis.rotationMetadata.normalizationConvention,source.normalizationConvention);
                    testCase.verifyEqual(basis.rotationMetadata.eigenvalues,reshape(source.rotationEigenvalues(chosen,1),1,[]),AbsTol=0);
                    testCase.verifyEqual(basis.rotationMetadata.residuals.relativePencilResidual,source.rotationResiduals.relativePencilResidual(chosen,:),AbsTol=0);
                    testCase.verifySize(basis.rotationMetadata.residuals.relativePencilResidual,[numel(chosen) 1]);
                    sourceValues = source.F([-1;-.4;0]);
                    testCase.verifyEqual(basis.F([-1;-.4;0]),sourceValues(:,chosen)./scales,AbsTol=3e-12);
                    testCase.verifyEqual(basis.canonicalEndpointLabels,["surface","bottom"]);
                    testCase.verifyEmpty(basis.constructionRecord.selection.initialColumnIndices);
                    testCase.verifyFalse(ismember("endpointResponse",basis.capabilities.normalizationNames));
                    testCase.verifyError(@() basis.withNormalization("endpointResponse"),"IMModeBasis:UnavailableNormalization");
                    restored = basis.withNormalization("rotationRightMetric");
                    testCase.verifyEqual(restored.F([-1;-.4;0]),sourceValues(:,chosen),AbsTol=3e-12);
                end
            end
        end

        function sourcePhysicsAndHorizontalIdentityMustAgree(testCase)
            background = edgeBackground("freeSurface");
            source = background.analyticalSolution().geostrophicZeroAPVModesAtWavenumber(1.2);
            inferred = wrapSource(source,background,1:2);
            testCase.verifyEqual(inferred.horizontalKind,"fixedWavenumber");
            testCase.verifyEqual(inferred.k,source.k,AbsTol=0);
            info = sourceInfo(1:2);
            testCase.verifyError(@() IMModeBasis(source=source,background=background,family="geostrophicZeroAPV",modeInfo=info,normalization="depthRMSF",k=1.3),"IMModeBasis:IncompatibleWavenumber");
            wrongG = IMModeBackground.constantStratification(N0=1,zDomain=[-1 0],rotationRate=background.rotationRate,latitude=33,g=3,surfaceCondition="freeSurface");
            wrongF0 = IMModeBackground.constantStratification(N0=1,zDomain=[-1 0],rotationRate=1.1*background.rotationRate,latitude=33,g=2,surfaceCondition="freeSurface");
            wrongSurface = edgeBackground("rigidLid");
            for physical = {wrongG,wrongF0,wrongSurface}
                testCase.verifyError(@() wrapSource(source,physical{1},1:2),"IMModeBasis:IncompatibleBackground");
            end
            hydrostatic = IMModes.hydrostaticFModes(background,solver=IMSolverAnalytical(),nInternalModes=1);
            testCase.verifyError(@() hydrostatic.sourceCoordinateMap,"IMModeBasis:UnavailableCoordinateMap");
        end

        function emptySelectionsRetainTheStructuralSurfaceRequirement(testCase)
            background = edgeBackground("freeSurface");
            bottom = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=IMSolverAnalytical(),energyBoundary=IMGeneralizedEnergyBoundary(g0=Inf,gd=.4));
            absent = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,energyBoundary=IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf));
            for basis = {bottom,bottom.selectColumns([]),absent}
                current = basis{1};
                testCase.verifyFalse(ismember("rotateSurfaceBuoyancy",current.capabilities.rotationNames));
                testCase.verifyError(@() current.rotateSurfaceBuoyancy(),"IMModeBasis:SurfaceEndpointRequired");
            end
            surface = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=IMSolverAnalytical(),energyBoundary=IMGeneralizedEnergyBoundary(g0=.2,gd=Inf));
            emptySurface = surface.selectColumns([]);
            testCase.verifyTrue(ismember("rotateSurfaceBuoyancy",emptySurface.capabilities.rotationNames));
            rotated = emptySurface.rotateSurfaceBuoyancy();
            testCase.verifyEqual(rotated.nColumns,0);
            testCase.verifySize(rotated.sourceCoordinateMap,[1 0]);
            testCase.verifyEqual(rotated.canonicalEndpointLabels,"surface");
        end

        function aSchurTunedSignedNullDirectionCannotReceiveUnitEnergy(testCase)
            background = edgeBackground("freeSurface");
            canonical = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=IMSolverAnalytical(),energyBoundary=IMGeneralizedEnergyBoundary(g0=0,gd=.4));
            [H,B0,Bd] = physicalForms(canonical);
            K = H+.4*Bd;
            g0 = -(K(1,1)-K(1,2)^2/K(2,2))/B0(1,1);
            direction = [1;-K(1,2)/K(2,2)];
            signed = H+g0*B0+.4*Bd;
            majorant = H+abs(g0)*B0+.4*Bd;
            testCase.verifyLessThan(abs(direction'*signed*direction),1e-13*(direction'*majorant*direction));
            basis = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1.2,solver=IMSolverAnalytical(),energyBoundary=IMGeneralizedEnergyBoundary(g0=g0,gd=.4));
            rotated = basis.rotateBoundaryDepth();
            [~,index] = min(abs(rotated.rotationMetadata.eigenvalues));
            nullDirection = rotated.selectColumns(index);
            testCase.verifyError(@() nullDirection.withNormalization("generalizedEnergy"),"IMModeBasis:NullNormalization");
            testCase.verifyGreaterThan(norm(nullDirection.F([-1;-.4;0])),.01);
        end

        function aRigidNeutralNullStateCannotBeNormalizedFromPositiveRoundoff(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,surfaceCondition="rigidLid");
            k = .3;
            basis = IMModes.geostrophicZeroAPVModesAtWavenumber(background,k,solver=IMSolverSpectral(nEVP=8,coordinateKind="z"),energyBoundary=IMGeneralizedEnergyBoundary(g0=0,gd=0));
            % In a neutral rigid column F is constant and only the
            % difference of endpoint G values carries physical energy.
            expectedH = background.f0^4/(2*background.g^2*k^6*background.D)*[1 -1;-1 1];
            [integratedH,~,~] = physicalForms(basis);
            testCase.verifyEqual(integratedH,expectedH,RelTol=3e-13,AbsTol=3e-11);
            testCase.verifyEqual(basis.energyMatrix,expectedH,RelTol=3e-13,AbsTol=3e-11);
            rotated = basis.rotateBoundaryDepth();
            [~,index] = min(abs(rotated.rotationMetadata.eigenvalues));
            nullDirection = rotated.selectColumns(index);
            z = [-1;-.7;-.2;0];
            testCase.verifyLessThan(norm(nullDirection.F(z)),1e-12);
            testCase.verifyGreaterThan(norm(nullDirection.G(z)),1);
            testCase.verifyLessThan(max(nullDirection.G(z))-min(nullDirection.G(z)),1e-12);
            testCase.verifyError(@() nullDirection.withNormalization("generalizedEnergy"),"IMModeBasis:NullNormalization");
        end
    end
end

function background = edgeBackground(surface)
background = IMModeBackground.constantStratification(N0=1,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,surfaceCondition=surface);
end

function info = sourceInfo(columns)
info = table("provided:"+string(columns(:)),repmat("endpointResponse",numel(columns),1),columns(:),VariableNames=["columnLabel","branch","branchIndex"]);
end

function basis = wrapSource(source,background,columns)
basis = IMModeBasis(source=source,background=background,family="geostrophicZeroAPV",modeInfo=sourceInfo(columns),columns=columns,normalization="depthRMSF",scientificInput=struct(energyBoundary=IMGeneralizedEnergyBoundary(g0=-3,gd=.4)));
end

function [H,B0,Bd] = physicalForms(basis)
[z,w] = legpts(64,basis.zDomain); w = w(:);
F = basis.F(z); G = basis.G(z);
f = basis.background.f0; g = basis.background.g; k = basis.k;
Fs = basis.F(basis.zDomain(2));
H = (k^2*(F'*(w.*F))+(f^2/g^2)*(G'*(w.*basis.background.N2(z).*G)))/(2*k^4);
if basis.background.surfaceCondition == "freeSurface", H = H+f^2/(2*g*k^4)*(Fs'*Fs); end
response = basis.G(basis.zDomain([2 1]).');
if basis.background.surfaceCondition == "freeSurface", response(1,:) = response(1,:)-Fs; end
B0 = f^2/(2*g^2*k^4)*(response(1,:)'*response(1,:));
Bd = f^2/(2*g^2*k^4)*(response(2,:)'*response(2,:));
end
