classdef IMModeBasisProvenanceTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            originalPath = path;
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            testCase.addTeardown(@() path(originalPath));
        end
    end

    methods (Test)
        function generalizedEnergyBoundaryRemainsIndependentOfPhysicalSurface(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            energy = IMGeneralizedEnergyBoundary(g0=-.04,gd=.03);
            basis = IMModes.geostrophicAPVModes(background,energyBoundary=energy,nPositiveModes=2,solver=IMSolverSpectral(nEVP=64));
            record = basis.scientificSpecification;
            testCase.verifyEqual(record.background.surfaceCondition,"rigidLid")
            testCase.verifyEqual(record.background.f0,2*background.rotationRate*sind(background.latitude))
            testCase.verifyEqual(record.requested.energyBoundary,struct(g0=-.04,gd=.03))
            testCase.verifyEqual(record.resolved.energyBoundary.g0,-.04)
            testCase.verifyEqual(record.resolved.energyBoundary.gd,.03)
            testCase.verifyEqual(record.resolved.evp.parameters.surfaceBoundary,"rigidLid")
            testCase.verifyEqual(record.resolved.evp.formulation,"F")
            testCase.verifyEqual(record.resolved.horizontal.kind,"independent")
            testCase.verifyEmpty(record.resolved.horizontal.k)
        end

        function resolvedDefaultAlphaIsDistinctFromRequestedDefault(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            basis = IMModes.geostrophicGeneralizedPotentialEnstrophyModes(background,1e-3,nModes=2,solver=IMSolverSpectral(nEVP=64));
            record = basis.scientificSpecification;
            testCase.verifyEmpty(record.requested.alpha0)
            testCase.verifyEqual(record.resolved.evp.parameters.alpha0Source,"stratification")
            testCase.verifyEqual(record.resolved.evp.parameters.bEffective,250,AbsTol=5e-12)
            testCase.verifyEqual(record.resolved.evp.parameters.alpha0,background.f0^2/250,RelTol=2e-14)
            testCase.verifyEqual(record.resolved.evp.parameters.alphaD,Inf)
            testCase.verifyEqual(record.resolved.horizontal.k,1e-3)
            testCase.verifyTrue(record.resolved.pairings.F.hasInnerProduct)
        end

        function constructionUsesConfiguredCoordinateWithoutClaimingConvergence(testCase)
            background = IMModeBackground.exponentialStratification(N0=.01,b=500,zDomain=[-1000 0],surfaceCondition="rigidLid");
            for coordinate = ["z","wkb","bandwidth"]
                basis = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=2,solver=IMSolverSpectral(nEVP=48,coordinateKind=coordinate));
                record = basis.constructionRecord;
                testCase.verifyEqual(record.requested.solver.coordinateKind,coordinate)
                testCase.verifyEqual(record.resolved.solver.coordinateKind,coordinate)
                testCase.verifyEqual(record.resolved.solver.nEVP,48)
                testCase.verifySize(record.resolved.solver.zNative,[48 1])
                testCase.verifyTrue(record.resolved.solver.nativeRepresentationAvailable)
                testCase.verifyEqual(record.selection.request.nPositiveModes,2)
                testCase.verifyNumElements(record.selection.sourceEigenpairIndices,2)
                testCase.verifyEqual(record.evidence.solveAccuracy.status,"unverified")
                testCase.verifyEqual(record.evidence.referenceConvergence.status,"unverified")
                testCase.verifyFalse(isfield(record,"recovery"))
            end
        end

        function finiteDifferenceRequestedGridIsReproducible(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            z = linspace(-1000,0,81).';
            basis = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=2,solver=IMSolverFiniteDifference(z=z));
            requested = basis.constructionRecord.requested.solver;
            actual = basis.constructionRecord.resolved.solver;
            testCase.verifyEqual(requested.className,"IMSolverFiniteDifference")
            testCase.verifyEqual(actual.className,"IMSolverFiniteDifference")
            testCase.verifyEqual(actual.coordinateKind,"z")
            testCase.verifyEqual(sort(requested.zNative),z)
            testCase.verifyEqual(actual.zNative,requested.zNative)
        end

        function emptyEVPKeepsItsScientificAndConstructionRequests(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=0,solver=IMSolverSpectral(nEVP=48,coordinateKind="bandwidth"));
            testCase.verifyEqual(basis.nColumns,0)
            testCase.verifyEqual(basis.scientificSpecification.resolved.evp.name,"waveModesAtWavenumber")
            testCase.verifyEqual(basis.scientificSpecification.requested.k,1e-3)
            testCase.verifyEmpty(basis.scientificSpecification.resolved.columnLabels)
            testCase.verifyEqual(basis.constructionRecord.resolved.solver.nEVP,48)
            testCase.verifyEqual(basis.constructionRecord.resolved.solver.coordinateKind,"bandwidth")
            testCase.verifyFalse(basis.constructionRecord.resolved.solver.nativeRepresentationAvailable)
            testCase.verifyFalse(basis.constructionRecord.resolved.constructionSkipped)
            testCase.verifyEmpty(basis.constructionRecord.selection.sourceEigenpairIndices)
        end

        function analyticalAndBarotropicBypassesDescribeActualConstruction(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            waves = IMModes.waveModesAtWavenumber(background,1e-3,nInternalModes=2,solver=IMSolverAnalytical());
            testCase.verifyEqual(waves.constructionRecord.resolved.backendKind,"analytical")
            testCase.verifyEqual(waves.constructionRecord.resolved.sourceClass,"IMAnalyticalInternalModesBasis")
            testCase.verifyEqual(waves.constructionRecord.resolved.solver.coordinateKind,"analytical")
            testCase.verifyEmpty(waves.constructionRecord.resolved.solver.nEVP)
            barotropic = IMModes.hydrostaticFModes(background,nInternalModes=0,solver=IMSolverSpectral(nEVP=48,coordinateKind="bandwidth"));
            testCase.verifyEqual(barotropic.constructionRecord.resolved.backendKind,"analyticalBarotropic")
            testCase.verifyTrue(barotropic.constructionRecord.resolved.constructionSkipped)
            testCase.verifyEmpty(fieldnames(barotropic.constructionRecord.resolved.solver))
            testCase.verifyEqual(barotropic.scientificSpecification.resolved.profileBoundarySemantics,"analytical rigid-lid F0=1, G0=0 before normalization")
            testCase.verifyTrue(barotropic.scientificSpecification.resolved.includesAnalyticalBarotropic)
            withInternals = IMModes.hydrostaticFModes(background,nInternalModes=2,solver=IMSolverSpectral(nEVP=48));
            testCase.verifyEqual(withInternals.family,"hydrostaticF")
            testCase.verifyEqual(withInternals.scientificSpecification.resolved.evp.formulation,"G")
            testCase.verifyTrue(withInternals.scientificSpecification.resolved.includesAnalyticalBarotropic)
            internalsOnly = withInternals.selectColumns([3 2]);
            testCase.verifyFalse(internalsOnly.scientificSpecification.resolved.includesAnalyticalBarotropic)
        end

        function DDIntervalsPreserveFixedEndpointAndCoefficientRoles(testCase)
            background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
            intervals = [-1000 -800;-200 0];
            energy = IMGeneralizedEnergyBoundary(g0=.03,gd=Inf);
            labels = IMModes.neutralLabelModes(background,neutralIntervals=intervals,nModesPerInterval=[1 2],energyBoundary=energy,normalization="maxAbsG");
            record = labels.scientificSpecification.resolved;
            testCase.verifyEqual(record.intervals,intervals)
            testCase.verifyEqual(record.intervalBoundaryTypes,repmat("Dirichlet",2,2))
            testCase.verifyEqual(record.endpointDerivativeTrace,"insideDeclaredInterval")
            testCase.verifyEqual(record.energyBoundaryRole,"provenanceOnly")
            testCase.verifyEqual(record.energyBoundary.g0,.03)
            testCase.verifyEqual(record.energyBoundary.gd,Inf)
            testCase.verifyEqual(record.profileRelation,"APV=-f0*Gz")
            testCase.verifyFalse(record.pairings.G.physicalEnergy)
            testCase.verifyEqual(labels.G(intervals(:)),zeros(4,3))
            zeroCount = IMModes.neutralLabelModes(background,neutralIntervals=intervals,nModesPerInterval=0);
            testCase.verifyEqual(zeroCount.scientificSpecification.resolved.intervals,intervals)
            testCase.verifyEmpty(zeroCount.scientificSpecification.resolved.columnLabels)
        end

        function criticalIntervalsAndEndpointSelectionsRemainExplicit(testCase)
            rotationRate = 7.2921e-5;
            latitude = 33;
            f0 = 2*rotationRate*sind(latitude);
            background = IMModeBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0],rotationRate=rotationRate,latitude=latitude);
            kernel = IMModes.zeroDepthWaveModesAtWavenumber(background,1e-3,criticalIntervals=[-900 -700],nModesPerInterval=2);
            testCase.verifyEqual(kernel.scientificSpecification.resolved.intervalKind,"critical")
            testCase.verifyEqual(kernel.scientificSpecification.resolved.profileRelation,"F=D*Gz")
            testCase.verifyEqual(kernel.scientificSpecification.resolved.horizontal.k,1e-3)
            testCase.verifyEqual(kernel.constructionRecord.resolved.backendKind,"analyticalDDKernel")
            energy = IMGeneralizedEnergyBoundary(g0=.03,gd=Inf);
            endpoint = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1e-3,energyBoundary=energy,solver=IMSolverSpectral(nEVP=48));
            testCase.verifyEqual(endpoint.scientificSpecification.resolved.endpoints,"surface")
            testCase.verifyEqual(endpoint.scientificSpecification.resolved.energyBoundaryRole,"definesPairingAndEndpointCoordinates")
            empty = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1e-3,energyBoundary=IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf));
            testCase.verifyEmpty(empty.scientificSpecification.resolved.endpoints)
            testCase.verifyEqual(empty.scientificSpecification.resolved.energyBoundary.g0,Inf)
            testCase.verifyTrue(empty.constructionRecord.resolved.constructionSkipped)
            testCase.verifyEmpty(fieldnames(empty.constructionRecord.resolved.solver))
        end

        function selectedRecordsComposeWithoutLosingOriginalRequest(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            original = IMModes.waveModesAtFrequency(background,.003,nPropagatingModes=4,solver=IMSolverAnalytical());
            selected = original.selectColumns([4 2 1]).selectColumns([2 1]);
            testCase.verifyEqual(selected.scientificSpecification.requested,original.scientificSpecification.requested)
            testCase.verifyEqual(selected.constructionRecord.requested,original.constructionRecord.requested)
            testCase.verifyEqual(selected.constructionRecord.selection.initialColumnIndices,[2 4])
            testCase.verifyEqual(selected.scientificSpecification.resolved.columnLabels,original.columnLabels([2 4]))
            testCase.verifyEqual(selected.scientificSpecification.resolved.modeInfo,original.modeInfo([2 4],:))
            testCase.verifyEqual(selected.scientificSpecification.resolved.normalization.scales,original.scientificSpecification.resolved.normalization.scales([2 4]))
            testCase.verifyEqual(selected.scientificSpecification.resolved.horizontal.kSquared,original.kSquared([2 4]))
            testCase.verifyEqual(selected.constructionRecord.selection.activeSourceColumns,selected.scientificSpecification.resolved.sourceColumns)
            testCase.verifyEqual(selected.constructionRecord.selection.sourceEigenpairIndices,original.constructionRecord.selection.sourceEigenpairIndices)
            record = selected.scientificSpecification;
            record.resolved.columnLabels(:) = "callerChange";
            testCase.verifyEqual(selected.columnLabels,original.columnLabels([2 4]))
            testCase.verifyEqual(selected.scientificSpecification.resolved.columnLabels,original.columnLabels([2 4]))
        end
    end
end
