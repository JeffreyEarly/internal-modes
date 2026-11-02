function basis = rotateZeroAPVPencil(self,name,left,right,ordering,convention,g0,gd)
% Apply one qualified pencil map without recreating the continuous source.
if self.family ~= "geostrophicZeroAPV"
    error("IMModeBasis:UnavailableRotation","Zero-APV rotations require the geostrophicZeroAPV family.");
end
[change,eigenvalues,rightSigns,residuals] = IMGeostrophicZeroAPVFormTools.solveSymmetricPencil(left,right,self.nColumns,1,ordering,"IMModeBasis");
change = change(:,:,1);
basis = self;
basis.zeroAPVColumnMap_ = self.sourceCoordinateMap*change;
basis.zeroAPVRightGram_ = change'*right*change;
if any(~isfinite(basis.zeroAPVColumnMap_),"all") || any(~isfinite(basis.zeroAPVRightGram_),"all")
    error("IMModeBasis:UnrepresentableCoordinates","The pencil produces a nonfinite source map or right form. Choose representable coefficient coordinates.");
end
basis.zeroAPVRightGram_ = IMGeostrophicZeroAPVFormTools.symmetrize(basis.zeroAPVRightGram_);
basis.columns = zeros(1,0);
basis.scaleFactors_ = ones(1,self.nColumns);
basis.normSignatures = reshape(rightSigns,1,[]);
basis.normalization = "rotationRightMetric";
labels = name+":"+string(1:self.nColumns);
basis.columnLabels = labels;
basis.modeInfo = table(labels(:),repmat("rotation",self.nColumns,1),(1:self.nColumns).',VariableNames=["columnLabel","branch","branchIndex"]);
h0 = zeros(1,0);
if convention == "unitEndpointResponseMetric", h0 = 2*self.k^2*reshape(eigenvalues,1,[]); end
basis.rotationMetadata = struct(name=name,isRotated=true,eigenvalues=reshape(eigenvalues,1,[]),rightMetricSignatures=reshape(rightSigns,1,[]),leftMetricSignatures=sign(reshape(eigenvalues,1,[])).*reshape(rightSigns,1,[]),normalizationConvention=convention,residuals=residuals,h0=h0,parentColumnLabels=self.columnLabels,coordinateChange=change,leftMatrix=left,rightMatrix=right);
basis.constructionRecord.selection.initialColumnIndices = zeros(1,0);
if ~isempty(g0), basis.rotationMetadata.g0 = g0; end
if ~isempty(gd), basis.rotationMetadata.gd = gd; end
basis.scientificSpecification.resolved.normalization = struct(name=basis.normalization,scales=basis.scaleFactors_,normSignatures=basis.normSignatures);
basis = basis.refreshZeroAPVRecords();
basis.capabilities = basis.buildCapabilities(struct());
end
