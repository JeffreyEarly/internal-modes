---
layout: default
title: Discrete transforms
parent: Class documentation V2
nav_order: 4
has_children: true
permalink: /classes-v2/discrete-transforms
mathjax: true
---

Reference pages for exact point-limited scalar transforms, mode-root grids, fixed-point quadrature fitting, forward and back modal transforms, and retained-band Gram, leakage, and quadratic-aliasing assessment.

The basic workflow starts from a solved `IMBasisSet` and requests an exact number of physical sample points. `discreteTransform` searches the available mode-root grids, fits one quadrature-weight vector against the largest candidate modal band having exactly that point count, and assesses every leading prefix without refitting the weights. It returns the largest consecutive prefix passing the default normalized-Gram policy.

```matlab
[transform,assessment] = basisSet.discreteTransform(nPoints=32);
assessment.prefixDiagnostics
```

`IMDiscreteTransform` stores `forwardMatrix` and `inverseMatrix`; use `transformForward` to compute retained modal coefficients and `transformBack` to return those coefficients to the fixed sample points. `IMDiscreteTransformAssessment` stores the returned production transform, full candidate transform, nested `IMQuadratureWeightFit`, point and mode counts, physical mode labels, per-prefix diagnostics, and the policy that limited retention.

The Gram policy uses the normalized operator error

$$
\left\lVert S_N(\Gamma_N-\Gamma_{0,N})S_N\right\rVert_2,
\qquad
S_N=\operatorname{diag}\!\left(\left|\operatorname{diag}\Gamma_{0,N}\right|^{-1/2}\right).
$$

Its default tolerance is $$10^{-2}$$. Constant and representative exponential regression sweeps set this value: constant cosine rules remain near roundoff, ordinary exponential cases in physical, WKB, and density coordinates remain below about $$5\times10^{-3}$$, and the unresolved DCT-I Nyquist prefix has order-one error. Override `gramTolerance` when an application needs a different Parseval-error budget.

Enable rejected-mode leakage or scalar quadratic-aliasing policies by supplying positive tolerances. Their empty defaults disable them:

```matlab
[transform,assessment] = basisSet.discreteTransform(nPoints=32, ...
    leakageTolerance=1e-3,quadraticAliasingTolerance=1e-2);
```

Leakage checks rejected single modes through `nCheckModes`, which defaults to twice the candidate count. Quadratic aliasing compares sampled product projections with continuous projections integrated independently on the source solver grid. Both policies require a positive-definite target metric; signed targets retain Gram diagnostics but reject these norm-based options. Coupled hydrostatic `F/G` transforms and product channels are a separate milestone.

Explicit-grid workflows remain available. Supplying `nModes` makes the requested band strict: if any enabled policy rejects it, the method reports the limiting policy, first failing value, tolerance, and maximum accepted count instead of silently shortening the band.

```matlab
[transform,assessment] = basisSet.discreteTransform(z=z,nModes=8);
[transform,assessment] = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
```

Use `quadratureWeightsForPoints` directly when the fitted weights or a custom linear objective are the primary concern. This lower-level method returns `IMQuadratureWeightFit` directly and does not select a retained prefix:

```matlab
[weights,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
```

These transform diagnostics assess one sampled rule, not whether the continuous EVP itself was solved accurately. Spectral residuals, coefficient tails, boundary residuals, and refinement belong to the separate solver-quality workflow.
