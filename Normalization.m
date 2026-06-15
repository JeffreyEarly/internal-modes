classdef Normalization
    % Enumerate internal-mode normalization conventions.
    %
    % `Normalization` names standard scaling conventions installed by
    % internal-mode basis sets. Each enum value resolves to a basis-set
    % normalization rule; `IMInternalModesBasis.normalizationFactors`
    % evaluates the selected rule and returns a per-mode factor $$s_j$$.
    % Evaluated internal-mode variables divide both $$F_j$$ and $$G_j$$ by
    % that same factor:
    % $$V_j^{\mathrm{out}}(z)=V_j^{\mathrm{raw}}(z)/s_j.$$
    % Generic canonical basis sets use string rule names directly, while
    % internal-mode basis sets may use either strings or these enum values.
    %
    % The internal-mode values are:
    %
    % - `Normalization.unity`: unit norm under the active internal-mode inner product
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
