# Vertical sampling grids in V2

A model samples every mode family, at every horizontal wavenumber, on one fixed set of depths. `IMVerticalGrid` holds those depths and their quadrature weights. It is separate from the numerical coordinate of an eigenvalue solver: modes can be solved at much higher resolution and evaluated on the grid afterwards.

```matlab
background = IMModeBackground(N2=@(z) 1e-6+5e-5*exp(z/300),zDomain=[-4000 0],latitude=33);
grid = IMVerticalGrid.smoothedBandwidth(background,65);
grid.z          % increasing depths, both endpoints included
grid.weights    % nonnegative weights with sum(grid.weights) == background.D
grid.integrate(background.N2(grid.z))
```

## Grid rules

| Constructor | Points | Weights |
| --- | --- | --- |
| `smoothedBandwidth(background,nPoints)` | Clenshaw–Curtis points in the smoothed-bandwidth coordinate | Mapped Clenshaw–Curtis weights, scaled so a constant integrates exactly |
| `physicalDepth(zDomain,nPoints)` | Clenshaw–Curtis points in depth | Clenshaw–Curtis |
| `modeRoots(basisSet,nPoints)` | Endpoints and roots of a generating mode (Early, Lelong and Smith 2020) | Nonnegative full-depth fit to the leading modes; can be zero where every fitted mode vanishes |
| `IMVerticalGrid(z=...,weights=...)` | Supplied | Supplied |

Every grid records its rule in `grid.kind` and its parameters and diagnostics in `grid.design`.

## The smoothed-bandwidth default

The stretching uses the local squared-frequency bandwidth

$$B=\sqrt{(N^2+\tilde f^2-f_0^2)^2+4f_0^2\tilde f^2},\qquad \tilde f=2\Omega\cos\phi,$$

which stays positive in neutral layers when the column rotates. It is smoothed with a Helmholtz filter, $$Q-L_s^2Q_{zz}=B$$, with $$Q_z=0$$ at both endpoints. Points are Clenshaw–Curtis nodes in $$x=\int\sqrt Q\,dz$$, mapped back to depth. Smoothing changes only where the points go; the physical problem keeps its own stratification and rotation.

The default smoothing length is

$$L_s = 2.5\,\frac{D}{N_z-1}.$$

This is provisional. It comes from the joint modal-recovery study in the `wave-vortex-model-coordinate-grid` note, for a 1000 m column with a 200 m neutral surface layer. At 65 points the default recovers 151 physical plus 2 neutral-label amplitudes, compared with 127 on a piecewise grid split at the neutral-layer joins. Smoothing lengths below about 1.2 grid spacings fail to resolve the mandatory boundary modes. At 33 points the best global grids lie near 2.4–3.4 spacings, while a fixed length in metres does markedly worse, so the length scales with the grid spacing. Other stratifications, boundary choices and wavenumbers are still being validated.

Supply `breakpoints` at depths where `N2` is not smooth, for example the ends of a neutral layer, so the bandwidth is represented accurately. Supply `smoothingLength` to override the default.

## Limits of any fixed grid

For a given number of points, some modes cannot be resolved by any grid. The boundary-buoyancy (zero-APV) responses thin as the horizontal wavenumber grows, and for a 1000 m column with a neutral surface layer no 33-point grid, global or piecewise, resolves the surface response to 1%. Such cases call for more points; they should be reported rather than worked around by tuning the grid.
