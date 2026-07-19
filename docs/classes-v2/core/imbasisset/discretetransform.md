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
 [transform,weightFit] = discreteTransform(basisSet,options)
```
## Parameters
+ `options.z`  increasing physical sample points
+ `options.weights`  optional quadrature weights aligned with `z`
+ `options.nModes`  number of leading retained modes

## Returns
+ `transform`  scalar discrete Galerkin transform
+ `weightFit`  quadrature-fit diagnostics, or empty when weights are supplied

## Discussion

For retained normalized modes $$u_j$$ sampled at points $$z_i$$, this
method forms

$$
(A_{\mathrm i})_{ij}=\Phi_{ij}=u_j(z_i),\qquad
W_{\mathrm{int}}=\operatorname{diag}\!\left(r(z_i)w_i\right),
$$

then constructs the Galerkin forward matrix stored by
`IMDiscreteTransform`. Eigenvalue-dependent endpoint terms are included
when they depend only on a sampled endpoint value. Endpoint derivative
traces cannot be inferred from arbitrary point samples and are rejected.
When `weights` is omitted, `quadratureWeightsForPoints` chooses
nonnegative weights with exact full-depth coverage by normalized Gram
Frobenius fitting. This fitted path requires `lsqlin` from Optimization
Toolbox.

```matlab
[transform,weightFit] = basisSet.discreteTransform(z=z,nModes=8);
transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
coefficients = transform.transformForward(values);
valuesFit = transform.transformBack(coefficients);
```

When `weights` is omitted, the optional second output preserves the fit
diagnostics and geometric comparison used to build the transform:

```matlab
[weightFit.residualNorm weightFit.geometricResidualNorm]
[weightFit.transform.relativeGramOperatorError weightFit.geometricTransform.relativeGramOperatorError]
```

In this case `transform` is the same fitted transform stored in
`weightFit.transform`. Supplying `weights` bypasses fitting, so the
optional `weightFit` output is empty.
