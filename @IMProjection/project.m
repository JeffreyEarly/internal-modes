function coefficients = project(self,values)
% Apply prepared Galerkin recovery or the supplied physical source system.
%
% Galerkin recovery and the continuous physical dual are distinct. Use
% projectPhysical to apply the continuous-normalization functional without
% its sampled-overlap correction. A failed recovery has no fallback inverse.
%
% - Topic: Apply projections
% - Declaration: coefficients = project(self,values)
% - Parameter values: finite real or complex sample-by-profile matrix
% - Returns coefficients: column-count-by-profile modal coefficients
arguments (Input)
    self (1,1) IMProjection
    values (:,:) double {mustBeFinite}
end
arguments (Output)
    coefficients (:,:) double
end
if size(values,1) ~= self.sampleCount
    error("IMProjection:InvalidShape","values must have one row per sample.");
end
if self.projectionKind == "prescribedDual"
    if self.physicalDualStatus ~= "available"
        error("IMProjection:UnavailablePhysicalDual","The supplied physical coefficient system is unavailable (%s); raw pairings remain available.",self.physicalDualStatus);
    end
elseif isempty(self.recovery) || ~self.recovery.success
    error("IMProjection:UnavailableRecovery","Sample recovery is unprepared or unresolved. Inspect recovery diagnostics, or use projectPhysical for the separately available continuous dual.");
end
coefficients = self.forwardMatrix*values;
end
