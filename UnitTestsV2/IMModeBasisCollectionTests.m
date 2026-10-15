classdef IMModeBasisCollectionTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function wavePagesPreserveTheStandardConventionAndScientificIdentity(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            first = IMModes.waveModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),nInternalModes=3);
            second = IMModes.waveModesAtWavenumber(background,2e-3,solver=IMSolverAnalytical(),nInternalModes=3);
            collection = IMBasisCollection({first,second},basisIndex=[2 1 2]);
            testCase.verifyEqual(collection.kappa,[2e-3 1e-3 2e-3]);
            testCase.verifyEqual(collection.metadata(1).family,first.family);
            testCase.verifyEqual(collection.metadata(1).columnLabels,first.columnLabels);
            testCase.verifyEqual(collection.metadata(1).normalization,"waveEnergy");
            testCase.verifyEqual(collection.metadata(1).variables,["F", "G"]);
            testCase.verifyEqual(collection.metadata(1).derivativeOrders,{0,0:2});
            z = linspace(-1000,0,17).';
            pages = [3 2 3 1];
            columns = [3 1 3];
            F = collection.evaluate(z,variable="F",pages=pages,columns=columns,sampleChunkSize=4,pageChunkSize=2);
            Gz = collection.evaluate(z,variable="G",derivativeOrder=1,pages=pages,columns=columns,sampleChunkSize=3,pageChunkSize=1);
            testCase.verifyEqual(F,background.D*Gz,AbsTol=2e-14);
            for j = 1:numel(pages)
                basis = collection.bases{collection.basisIndex(pages(j))};
                expected = basis.F(z);
                testCase.verifyEqual(F(:,:,j),expected(:,columns),AbsTol=0);
            end
            testCase.verifyError(@() collection.evaluate(z,variable="F",derivativeOrder=1),"IMBasisCollection:UnsupportedDerivative");
        end

        function streamedFacadePagesKeepRepeatedPositionsAndBoundedChunks(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            basis = IMModes.waveModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),nInternalModes=2);
            collection = IMBasisCollection({basis},basisIndex=[1 1 1]);
            z = linspace(-1000,0,19).';
            pages = [3 1 3 2 1];
            columns = [3 1];
            output = zeros(numel(z),numel(columns),numel(pages));
            coverage = zeros(numel(z),numel(pages));
            callCount = 0;
            collection.evaluateChunks(z,@consume,variable="G",derivativeOrder=2,pages=pages,columns=columns,sampleChunkSize=4,pageChunkSize=2);
            expected = basis.Gzz(z);
            testCase.verifyEqual(output,repmat(expected(:,columns),1,1,numel(pages)),AbsTol=0);
            testCase.verifyEqual(coverage,ones(size(coverage)));
            testCase.verifyEqual(callCount,15);

            function consume(values,rows,positions,pageIDs)
                testCase.verifyLessThanOrEqual(numel(rows),4);
                testCase.verifyLessThanOrEqual(numel(positions),2);
                testCase.verifyEqual(pageIDs,pages(positions));
                output(rows,:,positions) = values;
                coverage(rows,positions) = coverage(rows,positions)+1;
                callCount = callCount+1;
            end
        end

        function endpointFacadesUseScalarEvaluatorsAcrossCollectionPages(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            energy = IMGeneralizedEnergyBoundary(g0=0,gd=0);
            basis = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),energyBoundary=energy);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            testCase.verifyEqual(collection.metadata.columnKind,"endpoint");
            testCase.verifyEqual(collection.metadata.endpointLabels,["surface", "bottom"]);
            z = [-1000;-700;-100;0];
            for variable = ["F", "G"]
                expected = basis.evaluate(z,variable);
                sampled = collection.evaluate(z,variable=variable,pages=[2 1 2],columns=[2 1],sampleChunkSize=2,pageChunkSize=2);
                testCase.verifyEqual(sampled,repmat(expected(:,[2 1]),1,1,3),AbsTol=0);
            end
            testCase.verifyError(@() IMBasisCollection({basis},sourcePage=2),"IMBasisCollection:InvalidSourcePage");
        end

        function emptyEndpointInventoryRemainsAnEmptyContinuousBasis(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            energy = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf);
            basis = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1e-3,energyBoundary=energy);
            collection = IMBasisCollection({basis},basisIndex=[1 1]);
            testCase.verifyEqual(basis.nColumns,0);
            testCase.verifyEmpty(collection.metadata.columnLabels);
            testCase.verifySize(collection.evaluate([-1000;-500;0],variable="F"),[3 0 2]);
        end

        function neutralAPVProjectionPreservesTheFacadeNormalization(testCase)
            background = IMModeBackground(N2=@(z) 1e-4*double(z < -200),zDomain=[-1000 0]);
            energy = IMGeneralizedEnergyBoundary(g0=0,gd=Inf);
            depth = IMModes.neutralLabelModes(background,neutralIntervals=[-200 0],energyBoundary=energy,nModesPerInterval=3);
            maximum = IMModes.neutralLabelModes(background,neutralIntervals=[-200 0],energyBoundary=energy,nModesPerInterval=3,normalization="maxAbsG");
            collection = IMBasisCollection({depth,maximum},basisIndex=[2 1 2]);
            testCase.verifyEqual(collection.metadata(1).family,"neutralLabel");
            testCase.verifyEqual(collection.metadata(1).variables,["G", "APV"]);
            testCase.verifyEqual(collection.metadata(1).derivativeOrders,{0:2,0});
            z = linspace(-200,0,401).';
            weights = .5*ones(size(z));
            weights([1 end]) = .25;
            coefficients = [2;-.5;1.25];
            for page = [1 2]
                basis = collection.bases{collection.basisIndex(page)};
                projection = collection.projectionOnGrid(z,weights,page=page);
                explicit = basis.projectionOnGrid(z,weights);
                testCase.verifyEqual(projection.forwardMatrix,explicit.forwardMatrix,AbsTol=0);
                testCase.verifyEqual(projection.columnLabels,basis.columnLabels);
                testCase.verifyEqual(projection.forwardMatrix*basis.APV(z)*coefficients,coefficients,AbsTol=2e-13);
                selected = collection.projectionOnGrid(z,weights,page=page,variable="APV",columns=[3 1]);
                samples = basis.APV(z);
                testCase.verifyEqual(selected.forwardMatrix*samples(:,[3 1]),eye(2),AbsTol=2e-13);
            end
            sampled = collection.evaluate(z,variable="APV",pages=[3 1 2],columns=[3 1]);
            testCase.verifyEqual(sampled(:,:,1),sampled(:,:,2),AbsTol=0);
            testCase.verifyError(@() collection.projectionOnGrid(z,weights,variable="G"),"IMModeBasis:UnavailableProjection");
        end

        function collectionDoesNotInventAProjectionForWaveProfiles(testCase)
            background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
            basis = IMModes.waveModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),nInternalModes=2);
            collection = IMBasisCollection({basis});
            z = linspace(-1000,0,11).';
            weights = 100*ones(size(z));
            weights([1 end]) = 50;
            testCase.verifyError(@() collection.projectionOnGrid(z,weights),"IMModeBasis:UnavailableProjection");
        end
    end
end
