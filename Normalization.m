classdef Normalization
    % Enumerate the supported modal normalization conventions.
    %
    % `Normalization` names the scaling conventions used by basis sets.
    % Each convention selects a rule in `evp.normalizations` that returns a
    % per-mode factor $$s_j$$. Evaluated modes divide every variable in mode
    % $$j$$ by that factor:
    % $$V_j^{\mathrm{out}}(z)=V_j^{\mathrm{raw}}(z)/s_j.$$
    %
    % The valid values are:
    %
    % - `Normalization.unity`: unit norm under the active EVP inner product
    % - `Normalization.kConstant`: unit `G` inner-product norm for fixed-$$K$$ wave modes
    % - `Normalization.omegaConstant`: unit `F` inner-product norm for fixed-$$\omega$$ wave modes
    % - `Normalization.uMax`: scale by $$\max_z |F_j(z)|$$
    % - `Normalization.wMax`: scale by $$\max_z |G_j(z)|$$
    % - `Normalization.surfacePressure`: scale by the raw surface value of $$F$$
    % - `Normalization.geostrophic`: hydrostatic geostrophic normalization
    %
    % ```matlab
    % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
    % solver = IMSolverSpectral(nEVP=128);
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet.normalization = Normalization.geostrophic;
    % G = basisSet.G(z);
    % ```
    %
    % - Topic: Configure normalization
    % - Declaration: classdef Normalization
    enumeration
        unity, kConstant, omegaConstant, uMax, wMax, surfacePressure, geostrophic
    end
end
