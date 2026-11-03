# Hydrostatic mode-root transform

`IMInternalModesDiscreteTransform.hydrostaticRoot` builds the vertical transform for the rigid-lid hydrostatic problem with non-active boundaries: a rigid lid and rigid bottom, with no independent boundary displacement. It returns the final sampling grid, aligned F/G operators, retained-band assessment and reusable continuous `IMHydrostaticModesBasis` in one call. Family problems and solvers create continuous bases; the transform factory coordinates the subsequent sampling and qualification steps.

```matlab
background = IMBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
[transform,grid,assessment,basis] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,33);
a = transform.transformForward(u,variable="F");   % coefficients of F0..Fr
u = transform.transformBack(a,variable="F");
assessment.search                                   % every fitted candidate band
```

## Modes

One G family is solved and reused:

$$-G_j''=\frac{N^2}{g h_j}G_j,\qquad G_j(z_b)=G_j(z_s)=0,\qquad F_j=h_jG_j',$$

using the common `depthRMSF` normalization, $$\int F_iF_j\,dz=D\delta_{ij}$$ and $$\int (N^2/g)G_iG_j\,dz=(D/h_j)\delta_{ij}$$. All fields of each internal column use the same scale; this changes coefficient units relative to lower-level legacy geostrophic normalization, while normalized Gram errors and retained-band decisions are invariant under consistent rescaling. The analytical barotropic sector $$F_0=1$$, $$G_0=0$$ is appended explicitly. The aligned family is $$F_0,\ldots,F_r$$ and $$G_0,\ldots,G_r$$: `transform.modeNumber` is `0:r`, `h(1)` is `Inf`, and $$G_0$$ is inactive, so its forward row is zero. This relation, $$F=hG'$$, differs from the wave family's $$F=DG'$$.

Before any fitting, the retained candidates, the omitted-check modes and the generating mode are compared with an independent solve at `refinementFactor` (default 1.5) times the resolution. Equivalent depths, F and G shapes and the root points must agree within `convergenceTolerance` (default 1e-6). Without a supplied solver the spectral resolution starts at max(128, 4 Nz) and doubles for at most five resolution levels, capped at 2048 unless the starting resolution is larger. At each level WKB is attempted first and physical z provides a fallback for strongly stretched backgrounds. Root-count failures are recorded as underresolved attempts rather than aborting the search. A supplied `IMSolverSpectral` is used as given with no coordinate fallback; failure reports the tested coordinate and resolution. The evidence is in `assessment.convergence`.

## Reusing continuous modes and choosing a grid

The returned `basis` is the same retained continuous family stored as `transform.basis`, with labels `barotropic:0`, `positive:1`, ..., `positive:r`. It evaluates F, G and derivatives at arbitrary physical depths. The larger solved family supplies both the root-generating mode and omitted-mode checks; it is distinct from the returned retained band. No second independent F EVP defines the profiles.

```matlab
F = basis.F(zEvaluation);
Gz = basis.Gz(zEvaluation);
reused = basis.discreteTransform(grid);  % fixed-grid Galerkin; no EVP or fit
selected = basis.selectColumns([1 3 5]); % preserve scales and identities
selectedTransform = selected.discreteTransform(grid);
```

The fixed-grid bridge currently supports rigid-lid hydrostatic F/G families only. The same per-variable Galerkin operators are available through `basis.projectionOnGrid(z,w,variable="F")` or `IMBasisCollection({basis}).projectionOnGrid(z,w,variable="G")`; all use the same continuous metric and labels. Other common families retain their existing lower-level transform paths. It always constructs Galerkin coefficient recovery; the root convenience factory additionally has the larger full-span family needed for truncated collocation.

To construct a continuous family and obtain candidate points independently:

```matlab
solver = IMSolverSpectral(nEVP=256,coordinateKind="z");
family = solver.solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=33); % F0 plus 32 internal columns
[zRoots,rootDesign] = family.modeRootPoints(32);
```

`modeRootPoints` requires the generating mode to be present and returns points, without another solve or a weight fit. `IMVerticalGrid` holds the final physical points, quadrature weights and construction provenance. The solver's numerical coordinate grid remains independent. For continuous construction, `nModes=1` retains only the barotropic column, while `nModes=0` returns an empty basis. The root factory's separate `nInternalModes=0` option retains a barotropic-only transform band; its G channel has no active modes.

## Points and weights

With $$N_z$$ points including both endpoints, the interior points are the roots of $$G_{N_z-1}$$ (Early, Lelong and Smith 2020). Roots always come from the G family, whichever form an equivalent problem might be solved in. The point set is fixed; the retained count $$r\le N_z-2$$ and the solver resolution are independent of it.

For each candidate band one nonnegative weight vector with $$\sum_i w_i=D$$ is fitted to the stacked normalized F and G Gram errors, with $$W_F=\mathrm{diag}(w)$$ and $$W_G=\mathrm{diag}(N^2(z_i)w_i/g)$$. The fit uses `lsqlin` from Optimization Toolbox. The returned `grid.weights` and `transform.weights` are the fit for the returned band; nothing is fitted twice. Zero weights are allowed and counted in the assessment; code that divides by weights must handle them.

## Acceptance

With $$\Gamma_V=\Phi_V^\mathsf{T}W_V\Phi_V$$, $$H_F=D I_{r+1}$$ and the active $$H_G=\mathrm{diag}(D/h_1,\ldots,D/h_r)$$, a band is accepted when both sampled Gram matrices have full rank and

$$\varepsilon_V=\left\|H_V^{-1/2}(\Gamma_V-H_V)H_V^{-1/2}\right\|_2\le\texttt{gramTolerance},\qquad V\in\{F,G\}.$$

$$\varepsilon_V$$ bounds the relative error in the energy of any combination of retained modes computed from samples. Exact coefficient round trips and acceptable conditioning do not imply small $$\varepsilon_V$$.

- **Strict:** `nInternalModes=r` fits that band and either accepts it or throws `IMInternalModesDiscreteTransform:HydrostaticRootBandRejected` with both errors and ranks. It never truncates.
- **Automatic:** bands are fitted from $$N_z-2$$ down to zero and the largest accepted band is returned. `assessment.limitingCriterion` names the check that rejected band $$r+1$$, or `"maximumRootGridBand"` when $$r=N_z-2$$.

Acceptance is linear. It does not qualify nonlinear products and is not a claim of an optimal grid or a maximum resolution.

## Operators

The transform separates two jobs. The fitted weights are the rule for energies and other integrals on the grid, and the energy errors above qualify the retained band. A forward operator, chosen with `forwardOperator`, turns samples into coefficients. Both operators recover retained modes exactly and are cached, so repeated transforms are dense multiplies.

| `forwardOperator` | Coefficients | Modes inside the root span but outside the band | Use |
| --- | --- | --- | --- |
| `"truncatedCollocation"` (default) | Rows $$0,\ldots,r$$ of the inverses of the square systems $$[F_0,\ldots,F_{N_z-1}]$$ on all points and $$[G_1,\ldots,G_{N_z-2}]$$ on the interior points | Produce no retained coefficients | Model transforms |
| `"galerkin"` | $$A_V=\Gamma_V^{-1}\Phi_V^\mathsf{T}W_V$$ by weighted QR of $$W_V^{1/2}\Phi_V$$ | Leak into retained coefficients | Weighted least-squares projection, for example of data |

Truncated collocation combines v4's separation of modes within the interpolation span with a retained band qualified for sampled energy accuracy. Retaining fewer modes than the span does not change which modes are separated exactly, only which coefficients are kept. The [shared recovery engine](SharedRecovery.md) prepares the scaled full interpolation span once using rank-revealing sample factorization, then selects the retained output rows. Its canonical policy is `extendedSpan`; the physical adapter name remains `truncatedCollocation`. Condition numbers remain in `assessment.collocation`. Modes at or beyond the span ($$F_{N_z}$$, $$G_{N_z-1}$$ and higher) alias under either operator, and `transform.forwardOperator(variable=V)` reports which one a transform uses. Channel `recovery` and `physicalProjection` accessors expose the distinct prepared numerical map and physical quadrature record.

Two related quantities are different from either operator and are not substituted for them:

- `modeProjectionFunctional` returns the raw pairings $$\Phi_V^\mathsf{T}W_Vx$$;
- dividing active F pairings by $$H_F$$ approximates $$H_F^{-1}\int F_j x\,dz$$; dividing active G pairings by $$H_G$$ approximates $$H_G^{-1}\int (N^2/g)G_j x\,dz$$. These physical continuous-dual projections differ from exact sampled recovery by the Gram error, even for modal input. The inactive barotropic G row remains zero.

## Diagnostics

Each row of `assessment.search` reports, for the selected forward operator, $$\varepsilon_F$$, $$\varepsilon_G$$, ranks and Gram condition numbers, the weight range and zero-weight count, the largest normalized modal mean $$|\sum_iw_iF_j(z_i)|/D$$, and the omitted-mode gains $$\|H_r^{1/2}A_V\Phi_{V,\mathrm{om}}H_\mathrm{om}^{-1/2}\|_2$$ for the next `nOmittedModes` (default 4) internal modes. Gains are diagnostics for that finite omitted set, not acceptance criteria. Near $$r=N_z-2$$ they are typically of order one: the root grid aliases the next modes fully onto the retained band.

## Constant stratification

For constant $$N$$ the roots are uniform, the fitted weights are the trapezoid rule, and the transform is the correctly scaled DCT-I for F and DST-I for G with $$r=N_z-2$$ and $$\varepsilon_V$$ at roundoff. That case admits an $$O(N_z\log N_z)$$ fast path; this implementation uses dense matrices throughout.

## Comparison with v4

WVM v4 used the same root points and square collocation span, then truncated its operators to the requested `Nj` (commonly a reduced band for antialiasing). It used moment weights $$w=F^{-\mathsf{T}}[D,0,\ldots,0]^\mathsf{T}$$, which are exact for modal means but not for products. On identical points and converged modes, the default truncated-collocation rows equal v4's. The differences are the qualified band and the fitted weights for energies. `ExamplesV2/HydrostaticRootTransformComparison.m` reproduces v4's algorithm and reports usable bands under each weight rule, energy errors, omitted-mode gains for both operators, nonmodal-input errors and timings. In the recorded comparison the fitted weights qualify an equal or larger band, truncated collocation shows no leakage from modes inside the span, and Galerkin leaks them at roughly 1–17% amplitude. Neither choice is assumed to win in every case.

## Scope

Rigid lid, rigid bottom and $$N^2>0$$ only. Active boundaries, free surfaces and neutral layers use the families and grids described in [vertical sampling grids](VerticalGrids.md). Wavenumber-dependent wave families on this grid, WVM adoption and nonlinear qualification are separate work.
