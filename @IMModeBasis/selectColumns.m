function basis = selectColumns(self,columns)
% Select existing continuous columns without solving or renormalizing.
%
% The requested order, scientific identities, normalization scales and
% spectra are preserved. This is useful for retaining a band after a larger
% family supplied the generating mode of a sampling grid.
%
% - Topic: Inspect modes
% - Parameter columns: distinct indices into this basis, in the desired order
% - Returns basis: immutable selected IMModeBasis
arguments
    self (1,1) IMModeBasis
    columns (1,:) double {mustBeInteger,mustBePositive}
end
if any(columns > self.nColumns) || numel(unique(columns)) ~= numel(columns)
    error("IMModeBasis:InvalidColumns","Select distinct existing basis columns.");
end
basis = self;
if self.family == "geostrophicZeroAPV"
    basis.zeroAPVColumnMap_ = self.zeroAPVColumnMap_(:,columns);
    basis.zeroAPVRightGram_ = self.zeroAPVRightGram_(columns,columns);
    rotation = self.rotationMetadata;
    rotation.coordinateChange = rotation.coordinateChange(:,columns);
    for name = ["eigenvalues","rightMetricSignatures","leftMetricSignatures","h0"]
        if ~isempty(rotation.(name)), rotation.(name) = rotation.(name)(columns); end
    end
    if isfield(rotation.residuals,"relativePencilResidual")
        rotation.residuals.relativePencilResidual = rotation.residuals.relativePencilResidual(columns,:);
    end
    basis.rotationMetadata = rotation;
else
    basis.columns = self.columns(columns);
end
basis.nColumns = numel(columns);
basis.modeInfo = self.modeInfo(columns,:);
basis.columnLabels = self.columnLabels(columns);
basis.scaleFactors_ = self.scaleFactors_(columns);
basis.normSignatures = self.normSignatures(columns);
if self.hasEquivalentDepths_, basis.equivalentDepths_ = self.equivalentDepths_(columns); end
if self.hasEigenvalues_, basis.eigenvalues_ = self.eigenvalues_(columns); end
if ~isempty(self.frequency), basis.frequency = self.frequency(columns); end
if self.horizontalKind == "perColumn", basis.kSquared = self.kSquared(columns); end
resolved = basis.scientificSpecification.resolved;
resolved.columnLabels = basis.columnLabels;
resolved.modeInfo = basis.modeInfo;
resolved.sourceColumns = basis.columns;
resolved.normalization.scales = basis.scaleFactors_;
resolved.normalization.normSignatures = basis.normSignatures;
resolved.horizontal = struct(kind=basis.horizontalKind,k=basis.k,kSquared=basis.kSquared);
resolved.includesAnalyticalBarotropic = basis.family == "hydrostaticF" && any(basis.modeInfo.branch == "barotropic");
basis.scientificSpecification.resolved = resolved;
selection = basis.constructionRecord.selection;
selection.activeSourceColumns = basis.columns;
selection.columnLabels = basis.columnLabels;
if ~isempty(self.constructionRecord.selection.initialColumnIndices)
    selection.initialColumnIndices = self.constructionRecord.selection.initialColumnIndices(columns);
end
basis.constructionRecord.selection = selection;
if self.family == "geostrophicZeroAPV", basis = basis.refreshZeroAPVRecords(); end
end
