function self = refreshZeroAPVRecords(self)
% Separate original construction directions from current normalized maps.
map = self.zeroAPVColumnMap_;
isOneToOne = all(sum(map ~= 0,1) == 1);
if isOneToOne
    [rows,~] = find(map ~= 0);
    self.columns = reshape(rows,1,[]);
    coordinateKind = "canonicalEndpointColumns";
else
    self.columns = zeros(1,0);
    coordinateKind = "linearCombinations";
end
resolved = self.scientificSpecification.resolved;
resolved.columnLabels = self.columnLabels;
resolved.modeInfo = self.modeInfo;
resolved.sourceColumns = self.columns;
resolved.canonicalEndpointLabels = self.canonicalEndpointLabels;
resolved.sourceCoordinateMap = self.sourceCoordinateMap;
resolved.columnCoordinateKind = coordinateKind;
resolved.rotationMetadata = self.rotationMetadata;
resolved.endpoints = self.canonicalEndpointLabels;
resolved.energyBoundary = self.zeroAPVEnergyBoundary_;
resolved.energyBoundary.activeEndpoints = self.canonicalEndpointLabels;
resolved.energyBoundaryRole = "definesPairingAndEndpointCoordinates";
resolved.pairings.endpointResponse = struct(kind="canonicalEndpointResponseForm",canonicalEndpointLabels=self.canonicalEndpointLabels,coordinateFrame="sourceCoordinateMap");
resolved.pairings.generalizedEnergy = struct(kind="signedZeroAPVEnergy",formula="H + g0 B0 + gd Bd",g0=self.zeroAPVEnergyBoundary_.g0,gd=self.zeroAPVEnergyBoundary_.gd,coordinateFrame="current normalized columns");
self.scientificSpecification.resolved = resolved;
selection = self.constructionRecord.selection;
selection.activeSourceColumns = self.columns;
selection.columnLabels = self.columnLabels;
selection.canonicalEndpointLabels = self.canonicalEndpointLabels;
selection.sourceDirections = self.zeroAPVColumnMap_;
selection.columnCoordinateKind = coordinateKind;
if ~isOneToOne || self.rotationMetadata.isRotated
    selection.initialColumnIndices = zeros(1,0);
end
self.constructionRecord.selection = selection;
end
