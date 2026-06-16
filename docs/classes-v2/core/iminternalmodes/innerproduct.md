---
layout: default
title: innerProduct
parent: IMInternalModes
grand_parent: Core
nav_order: 12
mathjax: true
---

#  innerProduct

Return the `F` or `G` inner-product recipe.


---

## Declaration
```matlab
 spec = innerProduct(evp,variable)
```
## Parameters
+ `variable`  optional variable name, `"F"` or `"G"`

## Returns
+ `spec`  struct with interior and endpoint inner-product terms

## Discussion

  For `G`, the interior weight is $$N^2/g$$. For `F`, the
  interior weight is one. The returned struct has fields
  `variable`, `interiorWeight`, `surfaceWeights`,
  `bottomWeights`, `endpointFunctionals`, `hasInnerProduct`,
  and `reason`. `hasInnerProduct` is true when the variable has
  a known standalone inner product. When it is false, Gram
  matrices, spectra, and inner-product normalization for that
  variable throw `IMInternalModesBasis:UnavailableInnerProduct`.
  Diagnostic variables use the value-only hydrostatic endpoint
  catalog only when `modeFamily` is `"geostrophic"` and a
  catalog row is known; other diagnostic inner products are
  unavailable until a family catalog is added.
