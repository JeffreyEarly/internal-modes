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
+ `spec`  struct with interior and endpoint metric terms

## Discussion

  For `G`, the interior weight is $$N^2/g$$. For `F`, the
  interior weight is one. The returned struct has fields
  `variable`, `interiorWeight`, `surfaceWeights`,
  `bottomWeights`, `status`, and `reason`. `status` is
  `"fixed"` or `"interiorOnly"` when a standalone Gram matrix is
  available. It is `"unknown"`, `"mixed"`, or
  `"eigenvalueDependent"` when the requested diagnostic
  variable does not yet have an installed fixed inner-product
  rule.
