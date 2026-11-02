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
basis.columns = self.columns(columns);
basis.nColumns = numel(columns);
basis.modeInfo = self.modeInfo(columns,:);
basis.columnLabels = self.columnLabels(columns);
basis.scaleFactors_ = self.scaleFactors_(columns);
basis.normSignatures = self.normSignatures(columns);
if self.hasEquivalentDepths_, basis.equivalentDepths_ = self.equivalentDepths_(columns); end
if self.hasEigenvalues_, basis.eigenvalues_ = self.eigenvalues_(columns); end
if ~isempty(self.frequency), basis.frequency = self.frequency(columns); end
end
