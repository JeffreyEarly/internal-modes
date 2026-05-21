classdef Normalization
    % Enumerate the supported modal normalization conventions.
    %
    % `Normalization` collects the normalization choices used throughout
    % the `internal-modes` hierarchy. Following Section 2.4 of Early,
    % Lelong, and Smith (2020), the main wave-mode norms are the
    % `unity` norm based on the active EVP inner product, the `kConstant`
    % norm based on the fixed-$$K$$ orthogonality relation, and the
    % `omegaConstant` norm based on the fixed-$$\omega$$ relation.
    %
    % The valid values are:
    %
    % - `Normalization.unity` for unit norm under the active EVP inner product
    % - `Normalization.kConstant` for the manuscript's $$K$$-constant norm
    % - `Normalization.omegaConstant` for the manuscript's
    %   $$\omega$$-constant norm
    % - `Normalization.uMax` to scale each mode so $$\max F_j = 1$$
    % - `Normalization.wMax` to scale each mode so $$\max G_j = 1$$
    % - `Normalization.surfacePressure` to scale by the raw surface value
    %   of $$F$$ so $$F_j^\mathrm{surfacePressure}(z_\mathrm{surface})=1$$
    % - `Normalization.geostrophic` for the near-geostrophic interior mode
    %   normalization used by some helper workflows
    %
    % ```matlab
    % im = InternalModes(rho, zIn, zOut, latitude);
    % im.normalization = Normalization.kConstant;
    % ```
    %
    % - Topic: Configure normalization and boundaries
    % - Declaration: classdef Normalization
    enumeration
        unity, kConstant, omegaConstant, uMax, wMax, surfacePressure, geostrophic
    end
end
