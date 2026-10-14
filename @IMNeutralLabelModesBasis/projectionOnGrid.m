function projection = projectionOnGrid(self,z,weights,options)
% Project sampled APV using the manuscript's continuous normalization.
%
% Accepts APV itself (or its source), not a displacement or a full state:
% $$A_j = [f_0^2 D m_j^2]^{-1}\sum_r w_r q_j(z_r)^* \mathrm{APV}(z_r).$$
% Supply quadrature weights for dz and resolve each neutral interval.
% Interval-local quadrature avoids integrating endpoint derivative jumps
% across adjacent stratified water. Weights are never fitted. Quadrature
% error remains in the coefficients, even for sampled modal superpositions.
%
% Each pairing is divided by its exact norm before constructing an
% IMProjection with an identity coefficient system, avoiding small-f0 rank
% thresholds. The physical APV Gram is the error majorant and is recorded
% in provenance; targetGramMatrix is the row-scaled identity system. As a
% prescribed dual, this object has no synthesis/Gram assessment capability.
% Rank and conditioning describe that identity system, not whether the
% supplied samples resolve APV. Check forwardMatrix*basis.APV(z) directly
% against the identity on the retained columns to assess modal recovery.
% Zero rotation is rejected because APV cannot identify any neutral label.
%
% - Topic: Project APV
% - Parameter z: strictly increasing physical sample column within zDomain
% - Parameter weights: fixed quadrature weights for dz, one per sample
% - Parameter options.columns: distinct retained column indices, all by default
% - Returns projection: prescribed APV dual mapping samples to length amplitudes
arguments (Input)
    self IMNeutralLabelModesBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite}
    options.columns (1,:) double {mustBeInteger,mustBePositive} = []
end
arguments (Output)
    projection (1,1) IMProjection
end
if self.f0 == 0
    error("IMNeutralLabelModesBasis:ZeroRotationProjection","APV is identically zero at f0=0 and cannot identify neutral-label amplitudes.");
end
if numel(z) < 2 || any(diff(z) <= 0)
    error("IMNeutralLabelModesBasis:InvalidGrid","z must contain at least two strictly increasing points.");
end
if numel(weights) ~= numel(z) || ~any(weights ~= 0)
    error("IMNeutralLabelModesBasis:InvalidWeights","Supply one weight per point and at least one nonzero weight.");
end
columns = options.columns;
if isempty(columns)
    columns = 1:numel(self.m_n);
elseif any(columns > numel(self.m_n)) || numel(unique(columns)) ~= numel(columns)
    error("IMProjection:InvalidColumns","Select distinct existing neutral-label columns.");
end
% Cancel one factor of f0 analytically rather than forming 1/f0^2.
pairing = (-self.Gz(z).'./(self.f0*diff(self.zDomain)*self.m_n(:).^2)).*weights.';
if any(~isfinite(pairing(:)))
    error("IMNeutralLabelModesBasis:UnrepresentableProjection","The APV dual exceeds floating-point range for this rotation and geometry.");
end
recipe = self.projectionRecipe();
majorant = recipe.majorantGramMatrix(columns,columns);
if any(diag(majorant) <= 0) || any(~isfinite(majorant(:)))
    error("IMNeutralLabelModesBasis:UnrepresentableProjection","The physical APV norm is outside floating-point range for this rotation and geometry.");
end
projection = IMProjection.fromPrescribedDual(pairing(columns,:),eye(numel(columns)),eye(numel(columns)),majorantGramMatrix=majorant,columnLabels=self.columnLabels(columns),provenance=recipe.provenance);
end
