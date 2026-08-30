---
layout: default
title: discreteTransform
parent: IMBasisSet
grand_parent: Core
nav_order: 4
mathjax: true
---

#  discreteTransform

Build and assess a scalar Galerkin transform on a fixed quadrature rule.


---

## Declaration
```matlab
 [transform,assessment] = discreteTransform(basisSet,options)
```
## Parameters
+ `options.nPoints`  exact requested point count for a mode-root grid
+ `options.z`  increasing explicit physical sample points
+ `options.weights`  optional fixed quadrature weights aligned with `z`
+ `options.nModes`  optional strict number of leading candidate modes for explicit `z`
+ `options.gramTolerance`  nonnegative normalized-Gram operator-error tolerance
+ `options.leakageTolerance`  optional positive rejected-mode leakage tolerance
+ `options.quadraticAliasingTolerance`  optional positive scalar quadratic-aliasing tolerance
+ `options.nCheckModes`  optional number of source modes used by the leakage policy

## Returns
+ `transform`  retained production transform
+ `assessment`  fixed-rule candidate and retained-prefix diagnostics

## Discussion

A transform can be built from explicit physical points `z`, or from an
exact requested point count `nPoints`. The point-count workflow searches
the mode-root grids returned by `pointsFromModeRoots` and selects the
largest candidate modal prefix whose grid has exactly the requested
number of points. `nPoints` and `z` are mutually exclusive.

Unless `weights` are supplied with `z`, the quadrature weights are fitted
once for the complete candidate band by `quadratureWeightsForPoints`.
Every leading modal prefix is then assessed using those same physical
points and weights. The normalized-Gram policy is always active. Optional
rejected-mode leakage and scalar quadratic-aliasing policies are enabled
by supplying their positive tolerances. The returned production transform
contains the largest consecutive prefix accepted by every enabled policy.

If `nModes` is supplied with `z`, that band is strict: every requested
mode must pass every enabled policy, or construction throws an error.
With `nModes` omitted, the candidate band is the largest available prefix
that can be represented by the supplied points and may be reduced by the
policies. The `nPoints` workflow determines its candidate band from the
exact mode-root grid and therefore does not accept `nModes`.

For prefix $$N$$, define

$$
E_N=S_N(\Gamma_N-\Gamma_{0,N})S_N,\qquad
S_N=\operatorname{diag}\!\left(
\left|\operatorname{diag}\Gamma_{0,N}\right|^{-1/2}
\right).
$$

The always-enabled Gram policy accepts
$$\lVert E_N\rVert_2\leq\texttt{gramTolerance}$$. The shipped default
is $$10^{-2}$$. Regression sweeps over constant and representative
exponential stratifications in physical, WKB, and density coordinates
found constant cosine rules near roundoff and ordinary exponential rules
below approximately $$5\times10^{-3}$$; the factor-two margin preserves
those bands while still rejecting the order-one DCT-I Nyquist failure.

`leakageTolerance` enables rejected-mode leakage

$$
\ell_N=\max_{N<j\leq N_\mathrm{check}}
\frac{\lVert\Pi_N^\mathrm{discrete}u_j\rVert_\mu}
{\lVert u_j\rVert_\mu}.
$$

Its default check count is twice the candidate count. An explicit
`nCheckModes` must exceed the candidate band so every prefix has at least
one rejected comparison mode. `quadraticAliasingTolerance` enables

$$
q_N=\max_{1\leq i\leq j\leq N}
\frac{\left\lVert\Pi_N^\mathrm{discrete}(u_i u_j)-
\Pi_N^\mathrm{continuous}(u_i u_j)\right\rVert_\mu}
{\lVert u_i u_j\rVert_\mu}.
$$

The continuous quadratic projection is evaluated on the source solver's
inner-product grid, independently of the fixed transform rule. Leakage
and quadratic aliasing require a positive-definite target Gram matrix;
signed targets support Gram assessment only. `IMDiscreteTransformAssessment`
defines every table column and policy-result field.

```matlab
[transform,assessment] = basisSet.discreteTransform(nPoints=32);
[transform,assessment] = basisSet.discreteTransform(z=z,nModes=8);
[transform,assessment] = basisSet.discreteTransform(z=z,weights=w,nModes=8);
assessment.prefixDiagnostics
```
