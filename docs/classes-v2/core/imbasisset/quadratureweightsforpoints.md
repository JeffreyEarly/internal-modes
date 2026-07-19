---
layout: default
title: quadratureWeightsForPoints
parent: IMBasisSet
grand_parent: Core
nav_order: 21
mathjax: true
---

#  quadratureWeightsForPoints

Find quadrature weights for fixed physical sample points.


---

## Declaration
```matlab
 [weights,weightFit] = quadratureWeightsForPoints(basisSet,options)
```
## Parameters
+ `options.z`  increasing fixed physical sample points
+ `options.nModes`  number of leading retained modes
+ `options.objective`  `"normalizedGram"` or a custom least-squares callback
+ `options.nonnegative`  whether fitted weights must be nonnegative
+ `options.constrainDepth`  whether fitted weights must sum to the full depth

## Returns
+ `weights`  fitted quadrature weights aligned with `z`
+ `weightFit`  optional fitted-versus-geometric diagnostics

## Discussion

The points `z` are fixed. This method solves for one weight $$w_k$$ per
point so the sampled inner products of the first `nModes` retained modes
approximate their continuous Gram matrix. For target modal norms
$$C_i=(\Gamma_0)_{ii}$$ and fixed endpoint contribution
$$\Gamma_{\mathrm{endpoint}}$$, the default normalized least-squares
system is

$$
(A_{\mathrm{LS}})_{(i,j),k}
=\frac{r(z_k)\Phi_{ki}\Phi_{kj}}{\sqrt{|C_iC_j|}},
\qquad
(b_{\mathrm{LS}})_{(i,j)}
=\frac{(\Gamma_0-\Gamma_{\mathrm{endpoint}})_{ij}}
{\sqrt{|C_iC_j|}}.
$$

The fitted weights solve

$$
\min_w\left\|A_{\mathrm{LS}}w-b_{\mathrm{LS}}\right\|_2.
$$

By default they also satisfy

$$
w_k\geq0,
\qquad
\sum_k w_k=z_\mathrm{surface}-z_\mathrm{bottom}.
$$

Geometric control-volume weights provide the initial guess and reference
comparison. Unlike those geometric widths, fitted weights are algebraic
quadrature coefficients and need not remain positive or sum to the full
depth when the corresponding constraints are disabled. The optimizer
works with dimensionless weights $$x_k=w_k/D$$, where $$D$$ is the full
depth. Custom objectives continue to define $$Aw-b$$ in physical units.
This method requires `lsqlin` from Optimization Toolbox.

A custom objective is a function handle accepting a context struct and
returning a scalar struct with fields `A`, `b`, and optional `name`. The
context contains `z`, `modeNumber`, `normalization`, `inverseMatrix`,
`interiorWeight`, `targetGramMatrix`, `endpointGramMatrix`,
`geometricWeights`, `normalizedGramA`, and `normalizedGramB`.

```matlab
[weights,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
```
