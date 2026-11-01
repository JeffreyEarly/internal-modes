# Classical rigid-lid mode-root transform

`IMModes.classicalRootTransform` builds the vertical transform for the classical problem with non-active boundaries: a rigid lid and rigid bottom, with no independent boundary displacement. It returns the final sampling grid, the aligned F/G operators and an assessment of the retained band in one call.

```matlab
background = IMModeBackground(N2=@(z) 1e-4*exp(2*z/400),zDomain=[-1000 0],surfaceCondition="rigidLid");
[transform,grid,assessment] = IMModes.classicalRootTransform(background,33);
a = transform.transformForward(u,variable="F");   % coefficients of F0..Fr
u = transform.transformBack(a,variable="F");
assessment.search                                   % every fitted candidate band
```

## Modes

One G family is solved and reused:

$$-G_j''=\frac{N^2}{g h_j}G_j,\qquad G_j(z_b)=G_j(z_s)=0,\qquad F_j=h_jG_j',$$

with $$\int (N^2/g)G_iG_j\,dz=\delta_{ij}$$, so $$\int F_iF_j\,dz=h_j\delta_{ij}$$. The analytical barotropic sector $$F_0=1$$, $$G_0=0$$ is appended explicitly. The aligned family is $$F_0,\ldots,F_r$$ and $$G_0,\ldots,G_r$$: `transform.modeNumber` is `0:r`, `h(1)` is `Inf`, and $$G_0$$ is inactive, so its forward row is zero. This relation, $$F=hG'$$, differs from the wave creation layer's $$F=DG'$$.

Before any fitting, the retained candidates, the omitted-check modes and the generating mode are compared with an independent solve at `refinementFactor` (default 1.5) times the resolution. Equivalent depths, F and G shapes and the root points must agree within `convergenceTolerance` (default 1e-6). Without a supplied solver the WKB spectral resolution starts at max(128, 4 Nz) and doubles until the check passes; a supplied `IMSolverSpectral` is used as given. The evidence is in `assessment.convergence`.

## Points and weights

With $$N_z$$ points including both endpoints, the interior points are the roots of $$G_{N_z-1}$$ (Early, Lelong and Smith 2020). Roots always come from the G family, whichever form an equivalent problem might be solved in. The point set is fixed; the retained count $$r\le N_z-2$$ and the solver resolution are independent of it.

For each candidate band one nonnegative weight vector with $$\sum_i w_i=D$$ is fitted to the stacked normalized F and G Gram errors, with $$W_F=\mathrm{diag}(w)$$ and $$W_G=\mathrm{diag}(N^2(z_i)w_i/g)$$. The fit uses `lsqlin` from Optimization Toolbox. The returned `grid.weights` and `transform.weights` are the fit for the returned band; nothing is fitted twice. Zero weights are allowed and counted in the assessment; code that divides by weights must handle them.

## Acceptance

With $$\Gamma_V=\Phi_V^\mathsf{T}W_V\Phi_V$$, $$H_F=\mathrm{diag}(D,h_1,\ldots,h_r)$$ and the active $$H_G=I_r$$, a band is accepted when both sampled Gram matrices have full rank and

$$\varepsilon_V=\left\|H_V^{-1/2}(\Gamma_V-H_V)H_V^{-1/2}\right\|_2\le\texttt{gramTolerance},\qquad V\in\{F,G\}.$$

$$\varepsilon_V$$ bounds the relative error in the energy of any combination of retained modes computed from samples. Exact coefficient round trips and acceptable conditioning do not imply small $$\varepsilon_V$$.

- **Strict:** `nInternalModes=r` fits that band and either accepts it or throws `IMModes:ClassicalBandRejected` with both errors and ranks. It never truncates.
- **Automatic:** bands are fitted from $$N_z-2$$ downward and the largest accepted band is returned. `assessment.limitingCriterion` names the check that rejected band $$r+1$$, or `"maximumRootGridBand"` when $$r=N_z-2$$.

Acceptance is linear. It does not qualify nonlinear products and is not a claim of an optimal grid or a maximum resolution.

## Operators

`transformForward` applies the sampled Galerkin operator $$A_V=\Gamma_V^{-1}\Phi_V^\mathsf{T}W_V$$ on each active subspace, computed by weighted QR of $$W_V^{1/2}\Phi_V$$ and cached, so retained modes round-trip exactly and repeated transforms are dense multiplies. Two related quantities are different and are not substituted for it:

- `modeProjectionFunctional` returns the raw pairings $$\Phi_V^\mathsf{T}W_Vx$$;
- dividing those by $$H_V$$ approximates the prescribed continuous dual $$H_V^{-1}\int\Phi_Vx\,dz$$. It differs from Galerkin recovery by the Gram error, even for exactly modal input.

## Diagnostics

Each row of `assessment.search` reports $$\varepsilon_F$$, $$\varepsilon_G$$, ranks and Gram condition numbers, the weight range and zero-weight count, the largest normalized modal mean $$|\sum_iw_iF_j(z_i)|/\sqrt{Dh_j}$$, and the omitted-mode gains $$\|H_r^{1/2}A_V\Phi_{V,\mathrm{om}}H_\mathrm{om}^{-1/2}\|_2$$ for the next `nOmittedModes` (default 4) internal modes. Gains are diagnostics for that finite omitted set, not acceptance criteria. Near $$r=N_z-2$$ they are typically of order one: the root grid aliases the next modes fully onto the retained band.

## Constant stratification

For constant $$N$$ the roots are uniform, the fitted weights are the trapezoid rule, and the transform is the correctly scaled DCT-I for F and DST-I for G with $$r=N_z-2$$ and $$\varepsilon_V$$ at roundoff. That case admits an $$O(N_z\log N_z)$$ fast path; this implementation uses dense matrices throughout.

## Comparison with v4

WVM v4 used the same root points, square collocation with the auxiliary highest F mode, and moment weights $$w=F^{-\mathsf{T}}[D,0,\ldots,0]^\mathsf{T}$$, which are exact for modal means but not for products. `ExamplesV2/ClassicalRootTransformComparison.m` reproduces that algorithm from the same converged modes and reports usable bands, energy errors, omitted-mode gains, nonmodal-input errors and timings side by side. Neither method is assumed to win in every case. In the recorded comparison the fitted weights qualify an equal or larger band, while v4's square collocation, which inverts all $$N_z$$ modes before truncating, does not leak modes $$r+1,\ldots,N_z-1$$ into the retained coefficients; Galerkin projection onto a smaller band does, at roughly 1–17% amplitude in those profiles. Modes at or above $$N_z$$ alias under both.

## Scope

Rigid lid, rigid bottom and $$N^2>0$$ only. Active boundaries, free surfaces and neutral layers use the families and grids described in [vertical sampling grids](VerticalGrids.md). Wavenumber-dependent wave families on this grid, WVM adoption and nonlinear qualification are separate work.
