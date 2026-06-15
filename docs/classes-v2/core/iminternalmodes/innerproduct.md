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
  interior weight is one. Endpoint metric terms are attached to
  the solved formulation, because only that variable appears in
  the canonical endpoint condition. The returned struct has
  fields `variable`, `interiorWeight`, `surfaceWeights`, and
  `bottomWeights`.
