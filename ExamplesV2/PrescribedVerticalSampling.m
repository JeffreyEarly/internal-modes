function result = PrescribedVerticalSampling()
% Recover vertical coefficients from prescribed, staggered observations.
%
% F is observed at irregular interior points, while G is observed as true
% physical-depth cell averages. The locations need not be mode roots or
% include either endpoint. This example declares an unweighted fitting
% objective; it does not declare an observation-noise model.
arguments (Output)
    result (1,1) struct
end
background = IMModeBackground(N2=@(z) 1e-5 + 4e-5*exp(z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
basis = IMModes.hydrostaticFModes(background,solver=IMSolverSpectral(nEVP=64,coordinateKind="z"),nInternalModes=4);
zF = [-963;-812;-609;-377;-163;-29];
gCells = [-975 -821;-782 -565;-507 -302;-244 -36];
fSamples = IMVerticalSampling.pointValues("F",zF,observationLabels="F:"+string(1:numel(zF)));
gSamples = IMVerticalSampling.cellAverages("G",gCells,observationLabels="G:cell:"+string(1:size(gCells,1)));
sampling = IMVerticalSampling([fSamples.observations;gSamples.observations]);
[A,record] = sampling.evaluate(basis);

% The barotropic coefficient remains active: its F profile is constant,
% even though its G profile is exactly zero. Activity belongs to this
% shared coefficient problem, not to the G block alone.
transform = IMDiscreteTransform.fromObservationMatrices(A,columnLabels=record.columnLabels, ...
    activeColumnMask=true(1,basis.nColumns),fitFactor=eye(size(A,1)),metadata=record);
coefficients = reshape((1:30)+1i*(30:-1:1),basis.nColumns,2,3)/30;
observations = transform.transformBack(coefficients);
recovered = transform.transformForward(observations);
coefficientError = max(abs(recovered-coefficients),[],"all");

% A cell-center sample generally differs from the cell average.
midpointG = basis.G(mean(gCells,2));
cellAverageDifference = max(abs(A(numel(zF)+(1:size(gCells,1)),:)-midpointG),[],"all");
result = struct(basis=basis,sampling=sampling,transform=transform, ...
    coefficients=coefficients,observations=observations,recovered=recovered, ...
    coefficientError=coefficientError,cellAverageDifference=cellAverageDifference);
fprintf("%d prescribed observations, %d coefficients; round-trip error %.3g.\n",size(A,1),size(A,2),coefficientError);
fprintf("Cell means differ from cell-center samples by up to %.3g.\n",cellAverageDifference);
end
