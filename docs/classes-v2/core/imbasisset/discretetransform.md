---
layout: default
title: discreteTransform
parent: IMBasisSet
grand_parent: Core
nav_order: 4
mathjax: true
---

#  discreteTransform

Build a scalar Galerkin transform on fixed sample points.


---

## Declaration
```matlab
 transform = discreteTransform(basisSet,options)
```
## Parameters
+ `options.z`  increasing physical sample points
+ `options.increments`  optional quadrature increments aligned with `z`
+ `options.nModes`  number of leading retained modes

## Returns
+ `transform`  scalar discrete Galerkin transform

## Discussion

For retained normalized modes $$u_j$$ sampled at points $$z_i$$, this
method forms

$$
(A_{\mathrm i})_{ij}=\Phi_{ij}=u_j(z_i),\qquad
W_{\mathrm{int}}=\operatorname{diag}\!\left(r(z_i)\Delta z_i\right),
$$

then constructs the Galerkin forward matrix stored by
`IMDiscreteTransform`. Eigenvalue-dependent endpoint terms are included
when they depend only on a sampled endpoint value. Endpoint derivative
traces cannot be inferred from arbitrary point samples and are rejected.
When `increments` is omitted, `fitQuadrature` chooses nonnegative
increments with exact full-depth coverage by normalized Gram fitting.
This fitted path requires `lsqlin` from Optimization Toolbox.

```matlab
transform = basisSet.discreteTransform(z=z,nModes=8);
transform = basisSet.discreteTransform(z=z,increments=dz,nModes=8);
coefficients = transform.transformForward(values);
valuesFit = transform.transformBack(coefficients);
```
