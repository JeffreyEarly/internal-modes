function self = initializeZeroAPVCoordinates(self,scientificInput)
% Preserve canonical source rows separately from current profile columns.
if isempty(self.source)
    if self.nColumns ~= 0
        error("IMModeBasis:MissingCanonicalSource","A nonempty zero-APV basis requires a canonical numerical or analytical zero-APV source.");
    end
    endpoints = strings(1,0);
    if isfield(scientificInput,"endpoints"), endpoints = scientificInput.endpoints; end
    rotation = eye(numel(endpoints));
else
    if ~isa(self.source,"IMGeostrophicZeroAPVModesBasis") && ~isa(self.source,"IMAnalyticalGeostrophicZeroAPVModesBasis")
        error("IMModeBasis:MissingCanonicalSource","Use a numerical or analytical zero-APV source for this family.");
    end
    endpoints = reshape(self.source.endpoints,1,[]);
    rotation = self.source.rotationMatrix(:,:,1);
    if self.source.g ~= self.background.g || self.source.f0 ~= self.background.f0 || self.source.surfaceBoundary ~= self.background.surfaceCondition
        error("IMModeBasis:IncompatibleBackground","A zero-APV source and background must have matching g, f0 and physical surface convention. The developer wrapper also requires the same declared stratification.");
    end
    if self.horizontalKind == "unspecified" && isempty(self.k) && isempty(self.kSquared)
        self.horizontalKind = "fixedWavenumber";
        self.k = self.source.k;
    elseif self.horizontalKind ~= "fixedWavenumber" || ~isscalar(self.k) || self.k ~= self.source.k
        error("IMModeBasis:IncompatibleWavenumber","A zero-APV wrapper must preserve its source's fixed positive k. Omit an unassigned horizontal identity or supply the exact source wavenumber.");
    end
end
if self.horizontalKind ~= "fixedWavenumber" || ~isscalar(self.k) || self.k <= 0 || ~isempty(self.kSquared)
    error("IMModeBasis:InvalidHorizontalIdentity","A zero-APV basis requires its positive fixed horizontal wavenumber, including empty source sectors.");
end
self.canonicalEndpointLabels = endpoints;
self.zeroAPVSourceRotation_ = rotation;
self.zeroAPVColumnMap_ = rotation(:,self.columns);
self.zeroAPVRightGram_ = self.zeroAPVColumnMap_'*self.zeroAPVColumnMap_;
self.rotationMetadata.name = "boundaryNormalized";
self.rotationMetadata.normalizationConvention = "unitEndpointResponse";
self.rotationMetadata.parentColumnLabels = endpoints;
self.rotationMetadata.coordinateChange = self.zeroAPVColumnMap_;
boundary = struct(g0=Inf,gd=Inf);
if ismember("surface",endpoints), boundary.g0 = 0; end
if ismember("bottom",endpoints), boundary.gd = 0; end
if isfield(scientificInput,"energyBoundary")
    boundary.g0 = scientificInput.energyBoundary.g0;
    boundary.gd = scientificInput.energyBoundary.gd;
end
self.zeroAPVEnergyBoundary_ = boundary;
if ~isempty(self.source) && self.source.normalizationConvention ~= "unitEndpointResponse"
    self.rotationMetadata.name = self.source.rotationName;
    self.rotationMetadata.isRotated = true;
    self.rotationMetadata.normalizationConvention = self.source.normalizationConvention;
    self.rotationMetadata.eigenvalues = reshape(self.source.rotationEigenvalues(self.columns,1),1,[]);
    rightSigns = self.source.signatures(self.columns,1);
    if self.source.normalizationConvention == "unitEndpointResponseMetric", rightSigns = ones(size(rightSigns)); end
    self.rotationMetadata.rightMetricSignatures = reshape(rightSigns,1,[]);
    self.rotationMetadata.leftMetricSignatures = sign(self.rotationMetadata.eigenvalues).*reshape(rightSigns,1,[]);
    self.rotationMetadata.residuals = self.source.rotationResiduals;
    if isfield(self.rotationMetadata.residuals,"relativePencilResidual")
        self.rotationMetadata.residuals.relativePencilResidual = self.rotationMetadata.residuals.relativePencilResidual(self.columns,:);
    end
    if ~isempty(self.source.h0), self.rotationMetadata.h0 = reshape(self.source.h0(self.columns,1),1,[]); end
    self.zeroAPVRightGram_ = diag(rightSigns);
    if self.source.normalizationConvention == "unitEndpointResponseMetric"
        sourceRight = self.source.endpointResponseMetric;
        self.zeroAPVRightGram_ = sourceRight(self.columns,self.columns,1);
    end
end
end
