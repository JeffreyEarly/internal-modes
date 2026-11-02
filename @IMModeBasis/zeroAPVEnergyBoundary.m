function boundary = zeroAPVEnergyBoundary(self,g0,gd)
% Resolve signed coefficients without changing the canonical source sector.
if self.family ~= "geostrophicZeroAPV"
    error("IMModeBasis:UnavailableQuadraticForm","Generalized zero-APV energy requires the geostrophicZeroAPV family.");
end
boundary = self.zeroAPVEnergyBoundary_;
if ~isempty(g0), boundary.g0 = g0; end
if ~isempty(gd), boundary.gd = gd; end
for name = ["g0","gd"]
    value = boundary.(name);
    if ~isscalar(value) || ~isreal(value) || isnan(value) || value == -Inf
        error("IMModeBasis:InvalidEnergyCoefficient","g0 and gd must be finite signed scalars, or positive infinity for an absent endpoint.");
    end
end
if (isinf(boundary.g0) && ismember("surface",self.canonicalEndpointLabels)) || (isinf(boundary.gd) && ismember("bottom",self.canonicalEndpointLabels))
    error("IMModeBasis:InactiveEndpointRequired","An infinite coefficient is permitted only for a canonical endpoint absent from this source. Construct the desired source sector explicitly.");
end
end
