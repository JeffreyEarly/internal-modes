---
layout: default
title: innerProductNormFactor
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 10
mathjax: true
---

#  innerProductNormFactor

Return the `F` or `G` inner-product norm factor.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 factor = innerProductNormFactor(basisSet,iMode,options)
```
## Parameters
+ `iMode`  retained mode index
+ `options.variable`  `"F"` or `"G"`

## Returns
+ `factor`  raw inner-product norm factor

## Discussion

  This is the raw factor
  $$s_j=\sqrt{|\langle V_j,V_j\rangle|}$$ for `variable` equal
  to `F` or `G`. If `variable` is omitted, the solved
  formulation is used. The requested variable must have a known
  inner product. Custom normalization rules
  registered with `addNormalization` call this method.
  The factor is computed from raw, unnormalized modes before
  the active basis normalization is applied.
