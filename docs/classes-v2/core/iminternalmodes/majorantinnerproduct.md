---
layout: default
title: majorantInnerProduct
parent: IMInternalModes
grand_parent: Core
nav_order: 14
mathjax: true
---

#  majorantInnerProduct

Return the induced positive Hilbert-majorant recipe.


---

## Declaration
```matlab
 spec = majorantInnerProduct(evp,variable)
```
## Parameters
+ `variable`  optional variable name, `"F"` or `"G"`

## Returns
+ `spec`  positive interior and absolute-endpoint recipe

## Discussion

  The signed Pontryagin product returned by `innerProduct` is
  the physical pairing used for orthogonality and projection.
  Its natural $$L^2\oplus\mathbb C^s$$ coordinate decomposition
  induces a positive Hilbert product by retaining the positive
  interior weight and replacing every endpoint coefficient by
  its absolute value:

  $$
  (U,V)_+=\int w\,\overline{U}V\,dz+
  \sum_\ell |\alpha_\ell|\,
  \overline{L_\ell[U]}L_\ell[V].
  $$

  Use this recipe for magnitudes, error tolerances, and
  convergence diagnostics. Use `innerProduct` for signed
  invariants, projection functionals, and modal coefficients.
  The two recipes coincide when every endpoint coefficient is
  nonnegative.
