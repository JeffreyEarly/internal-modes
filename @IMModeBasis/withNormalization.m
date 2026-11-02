function basis = withNormalization(self,name)
% Change normalization using the existing continuous mode representation.
%
% The returned value preserves the scientific columns, spectra, physical
% background and original construction records. Every profile and its
% derivatives in a column use the same new positive divisor; signed metric
% signatures remain separate. Resolved normalization and, for zero-APV
% bases, the current source-coordinate map are updated. Rotation creation
% evidence and fixed construction directions retain their original meaning.
% Computing a new rule may evaluate the retained profiles and its fixed
% integration rule, but does not solve another EVP or fit sampling weights.
% Requesting the current normalization returns the existing value directly.
%
% ```matlab
% rmsModes = modes.withNormalization("depthRMSF");
% ```
%
% - Topic: Inspect modes
% - Declaration: basis = withNormalization(self,name)
% - Parameter name: one of the basis capability normalizationNames
% - Returns basis: immutable basis with the requested admissible normalization
arguments (Input)
    self (1,1) IMModeBasis
    name (1,1) string
end
arguments (Output)
    basis (1,1) IMModeBasis
end
basis = self;
if name == self.normalization
    return
end
basis.normalization = name;
[basis.scaleFactors_,basis.normSignatures] = basis.normalizationScales();
basis.scientificSpecification.resolved.normalization = struct(name=name,scales=basis.scaleFactors_,normSignatures=basis.normSignatures);
if self.family == "geostrophicZeroAPV", basis = basis.refreshZeroAPVRecords(); end
end
