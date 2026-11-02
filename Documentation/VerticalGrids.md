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

This is provisional. It comes from the joint modal-recovery study in the `wave-vortex-model-coordinate-grid` note. The historical counts hold only under that study's conditions: a 1000 m column with a 200 m neutral surface layer, latitude 45°, a 10 km wavelength, `g0 = 0` and `gd = -0.01`, data containing only the retained modes, and the study's joint overlap correction applied after the direct projectors, with fixed energy, enstrophy, boundary and noise tolerances. That study used the former generalized-energy-dependent neutral-label boundary policy, including a Neumann surface for the neutral cap. Its reported 65-point recovery of 151 physical plus two neutral-label amplitudes, against 127 on a piecewise grid split at the joins, is historical evidence rather than qualification of the current DD implementation. Smoothing lengths below about 1.2 grid spacings failed to resolve the mandatory boundary modes in those tests. `IMVerticalGrid` supplies only points and weights; ordinary manuscript projections on these points do not inherit that recovery result.

The current [neutral-label API](NeutralLabelModes.md) fixes both interval endpoint labels, independently of `g0` and `gd`. This DD correction preserves the provisional grid construction but selects no new retained-count limit. Portable joint-recovery qualification with the current shared engine belongs to [#58](https://github.com/JeffreyEarly/internal-modes/issues/58); local research requalification reports are separate from this provider delivery.

Retained-band recovery is a separate question from contamination by omitted modes. In the historical study, near the maximum retained counts, unresolved modes produced large false retained amplitudes on every tested grid, and retaining about 10% fewer modes per family reduced that aliasing sharply. This is not a current DD alias-rejection guarantee.

At 33 points, among the tested candidates and with the study's degree-five local reconstruction, the smallest boundary errors occur for smoothing lengths near 2.4–3.4 grid spacings, while a fixed length in metres does markedly worse. On that evidence the length is scaled with the grid spacing. Other stratifications, boundary choices and wavenumbers are still being validated.

Supply `breakpoints` at depths where `N2` is not smooth, for example the ends of a neutral layer, so the bandwidth is represented accurately. Supply `smoothingLength` to override the default.

## Limits of a fixed grid

A fixed number of points resolves only part of each family, and the boundary-buoyancy (zero-APV) responses thin as the horizontal wavenumber grows. In the 33-point comparison for the 1000 m column with a neutral surface layer, none of the tested grids (global smoothing lengths from 18 m to 2.4 km, physical depth, raw bandwidth and the piecewise benchmark) resolved the surface zero-APV response to the 1% shape tolerance under the study's degree-five reconstruction; the best reached 1.3% (piecewise) and 2.7% (smoothed). That finite sweep does not show that every 33-point grid or reconstruction must fail. More points are the natural next step in such a case, and the achieved boundary-mode accuracy should be reported rather than assumed.
