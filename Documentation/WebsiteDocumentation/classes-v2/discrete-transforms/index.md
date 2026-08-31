---
layout: default
title: Discrete transforms
parent: Class documentation V2
nav_order: 4
has_children: true
permalink: /classes-v2/discrete-transforms
mathjax: true
---

Reference pages for scalar and aligned internal-mode discrete transforms, mode-root grids, fixed-point shared quadrature fitting, and retained-band Gram, leakage, and quadratic-aliasing assessment.

## Scalar transforms

For a scalar `IMBasisSet`, `discreteTransform` searches the available mode-root grids for an exact requested point count, fits one quadrature rule against the largest matching candidate band, and assesses every leading prefix without refitting the weights:

```matlab
[transform,assessment] = basisSet.discreteTransform(nPoints=32);
coefficients = transform.transformForward(values);
valuesFit = transform.transformBack(coefficients);
```

`IMDiscreteTransform` stores the scalar forward and inverse matrices. `IMDiscreteTransformAssessment` stores the production prefix, candidate transform, fixed-rule prefix diagnostics, and active policies. `IMQuadratureWeightFit` preserves the full-band fitted rule and its geometric control-volume comparison.

The scalar normalized-Gram policy uses

$$
\left\lVert S_N(\Gamma_N-\Gamma_{0,N})S_N\right\rVert_2,
\qquad
S_N=\operatorname{diag}\!\left(\left|\operatorname{diag}\Gamma_{0,N}\right|^{-1/2}\right),
$$

with default tolerance $$10^{-2}$$. Optional rejected-mode leakage and scalar quadratic-aliasing policies require positive-definite targets. Signed scalar targets retain the Gram diagnostic only.

## Aligned internal-mode transforms

`IMInternalModesBasis` overrides both transform builders. One shared point-and-weight rule carries aligned $$F$$ and $$G$$ mode families, while each variable retains its own metric, active columns, target Gram matrix, forward projection, and diagnostics:

```matlab
[transform,assessment] = basisSet.discreteTransform( ...
    nPoints=24,variables=["F","G"]);

aG = transform.transformForward(GValues,variable="G");
FValues = transform.transformBack(aG,variable="F");
```

Omitting `variables` selects every directly representable channel in canonical order `F`, `G`. Explicit variables are validated and returned in that same order. A direct forward channel requires a known continuous inner product whose endpoint functionals are value-only, use the requested variable, and have every required endpoint present in `z`. Companion-variable and derivative endpoint traces remain available as diagnostics, but paired admissible-state inversion is outside this transform.

Let $$A_\mathrm{i}^{V}$$ be the sampled synthesis matrix for variable $$V$$, $$W_V$$ its sampled metric, and $$Q_V$$ the diagonal projector onto active family columns. A column is inactive only when its continuous norm and sampled variable column are both numerically zero, as for the barotropic $$G_0$$ column. The direct matrices obey

$$
\Gamma_V=(A_\mathrm{i}^{V})^\mathsf{T}W_VA_\mathrm{i}^{V},
\qquad
A_\mathrm{f}^{V}A_\mathrm{i}^{V}=Q_V.
$$

`modeProjectionFunctional(values,variable=V)` returns the unsolved pairings

$$
P_V(X)=(A_\mathrm{i}^{V})^\mathsf{T}W_VX,
$$

whereas `transformForward` applies the active Gram solve and returns coefficient arrays with inactive rows set to zero. `transformBack` always accepts a full family-sized coefficient array. `inverseMatrix` and `endpointValues` remain available for diagnostic-only variables even when no direct forward metric exists.

`IMInternalModesDiscreteTransform` also snapshots `modeNumber`, `h`, normalization, `zDomain`, depth, gravity, mode family, sampled `N2Values`, surface/bottom endpoint locations and traces, and basis metadata. It does not retain the source basis. In particular, geostrophic APV transforms preserve `g0`, `gd`, and `surfaceBoundary` in `problemMetadata`.

## Shared quadrature fitting

Use `quadratureWeightsForPoints` when fixed locations or custom fitting objectives are the primary concern:

```matlab
[weights,weightFit] = basisSet.quadratureWeightsForPoints( ...
    z=z,nModes=8,variables=["F","G"]);
```

For every requested channel, the builder compresses the normalized Gram Frobenius objective to its active upper-triangle mode pairs and stacks the systems without extra channel weighting:

$$
\min_w\sum_{V\in\mathcal V}
\left\lVert S_V\left(\Gamma_V(w)-\Gamma_{0,V}\right)S_V\right\rVert_\mathrm{F}^{2}.
$$

One nonnegative weight vector constrained to the full depth is fitted by default. `IMInternalModesQuadratureWeightFit` stores fitted and geometric specialized transforms, combined residual norms, per-variable residual norms, optimizer diagnostics, mode-pair indices, and row-to-variable provenance. Custom callbacks receive the stacked system and scalar-style contexts in `context.variableContexts.F` and `.G`.

## Common retained-band policies

The internal-mode Gram policy evaluates each requested variable independently and accepts only their common consecutive prefix. Leakage projects rejected $$F$$ modes into retained $$F$$ and rejected $$G$$ modes into retained $$G$$, then records the worst enabled variable. Coupled quadratic aliasing evaluates

$$
F_iF_j\rightarrow F,
\qquad F_iG_j\rightarrow G,
\qquad G_iG_j\rightarrow F.
$$

Products involving an identically zero source column are skipped. Continuous product projections use the solver's integration grid independently of the fixed transform points and weights. Leakage and quadratic-aliasing policies require positive-definite active target subspaces for every participating channel; signed channels use the Gram policy alone.

`IMInternalModesDiscreteTransformAssessment.prefixDiagnostics` records the worst variable or product channel at every prefix. `variablePrefixDiagnostics(variable=...)` exposes active counts, Gram and round-trip errors, rank, conditioning, and target definiteness for one channel. Supplying explicit `nModes` remains strict: a failing requested family band throws instead of silently shortening it.

These diagnostics assess the sampled transform rule. Continuous EVP residuals, coefficient tails, boundary residuals, and refinement remain part of the separate solver-quality workflow.

## Geostrophic APV and boundary composition

`IMGeostrophicTransform` consumes a generalized-energy APV `IMInternalModesDiscreteTransform` and a canonical `IMGeostrophicZeroAPVModesBasis`. Finite endpoint accelerations, including zero, select the active zero-APV coordinates; `Inf` omits that endpoint. The inputs must agree on stratification, domain, gravity, endpoint parameters, and free-surface or rigid-lid convention.

For APV eigendepth $$h_j$$ and positive horizontal wavenumber $$\kappa$$,

$$
\mu_\kappa^j=\kappa^2+\frac{f_0^2}{g h_j}.
$$

The transform rejects a retained mode when the relative separation of $$\mu_\kappa^j$$ from zero is at or below `muTolerance`. Its public `apvEndpointResponse` has dimensions `nEndpoints x nAPVModes x nK` and is assembled from the normalized APV endpoint traces:

$$
\mathbf r_q^{\kappa j}=-\frac{f_0}{g\mu_\kappa^j}
\begin{bmatrix}B_{\mathrm s}^j\\G_j(z_b)\end{bmatrix},
\qquad
B_{\mathrm s}^j=
\begin{cases}
G_j(0)-F_j(0), & \text{free surface},\\
G_j(0), & \text{rigid lid}.
\end{cases}
$$

Admissible-state arrays preserve explicit wavenumber pages and arbitrary trailing field dimensions:

| Quantity | Leading dimensions |
| --- | --- |
| sampled APV or sampled source | `nZ x nK x ...` |
| endpoint anomalies | `nEndpoints x nK x ...` |
| APV coefficients | `nAPVModes x nK x ...` |
| zero-APV coefficients | `nEndpoints x nK x ...` |

The state transform first projects APV through the `F` forward channel, then obtains boundary-normalized coefficients from

$$
\mathbf A_0^\kappa=-\frac{g\kappa^2}{f_0}
\left(\mathbf b^\kappa-\mathsf R_q^\kappa\mathbf A_q^\kappa\right).
$$

Its inverse reconstructs

$$
q=A_{\mathrm i}^{F}\mathbf A_q,
\qquad
\mathbf b^\kappa=\mathsf R_q^\kappa\mathbf A_q^\kappa
-\frac{f_0}{g\kappa^2}\mathbf A_0^\kappa.
$$

Generic source projection is deliberately separate. It uses `modeProjectionFunctional`, not `transformForward`, because the source equations require raw mode pairings before the APV Gram solve:

$$
S_q^j=\frac{1}{D}\int F_jS_\omega\,dz
-\frac{f_0}{D}\mathcal G_q^j[S_\eta].
$$

The canonical zero-APV source coefficients solve

$$
\widehat{\mathsf H}_g^\kappa\mathbf S_0^\kappa=\mathbf p_0^\kappa,
\qquad
\widehat{\mathsf H}_g^\kappa=\frac{2\kappa^2}{D}\mathsf H_g^\kappa,
$$

using the APV transform's points and weights together with the required endpoint source terms. Passing a rotated basis through `zeroAPVCoordinates` applies the inverse coordinate map on forward operations and the forward map on reconstruction, leaving the physical state and source projection unchanged.
