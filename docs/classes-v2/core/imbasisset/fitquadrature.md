---
layout: default
title: fitQuadrature
parent: IMBasisSet
grand_parent: Core
nav_order: 8
mathjax: true
---

#  fitQuadrature

Fit quadrature increments on fixed physical sample points.


---

## Declaration
```matlab
 fit = fitQuadrature(basisSet,options)
```
## Parameters
+ `options.z`  increasing fixed physical sample points
+ `options.nModes`  number of leading retained modes
+ `options.objective`  `"normalizedGram"` or a custom least-squares callback
+ `options.nonnegative`  whether fitted increments must be nonnegative
+ `options.constrainDepth`  whether fitted increments must sum to the full depth

## Returns
+ `fit`  fitted and geometric quadrature comparison

## Discussion

The default objective fits the normalized scalar Gram matrix. For target
norms $$C_i=(\Gamma_0)_{ii}$$ and a fixed endpoint contribution
$$\Gamma_{\mathrm{end}}$$, the ordered mode-pair residual system is

$$
A_{(i,j),k}=\frac{r(z_k)\Phi_{ki}\Phi_{kj}}{\sqrt{|C_iC_j|}},
\qquad
b_{(i,j)}=\frac{(\Gamma_0-\Gamma_{\mathrm{end}})_{ij}}
{\sqrt{|C_iC_j|}}.
$$

The fitted increments minimize $$\|A\Delta z-b\|_2$$. By default they
also satisfy $$\Delta z_k\geq0$$ and
$$\sum_k\Delta z_k=z_\mathrm{surface}-z_\mathrm{bottom}$$. The result
compares the fitted transform with geometric control-volume increments.
The optimizer solves for dimensionless increments
$$x_k=\Delta z_k/D$$, where $$D$$ is the full depth, while custom
objectives continue to define $$A\Delta z-b$$ in physical units.
This method requires `lsqlin` from Optimization Toolbox.

A custom objective is a function handle accepting a context struct and
returning a scalar struct with fields `A`, `b`, and optional `name`. The
context contains `z`, `modeNumber`, `normalization`, `basisMatrix`,
`interiorWeight`, `targetGramMatrix`, `endpointGramMatrix`,
`geometricIncrements`, `normalizedGramA`, and `normalizedGramB`.

```matlab
fit = basisSet.fitQuadrature(z=z,nModes=8);
transform = fit.fittedTransform;
```
