---
layout: default
title: majorantGramMatrix
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 15
mathjax: true
---

#  majorantGramMatrix

Return the positive Hilbert-majorant Gram matrix.


---

## Declaration
```matlab
 gram = majorantGramMatrix(basisSet,options)
```
## Parameters
+ `options.variable`  `"F"` or `"G"`
+ `options.zBounds`  integration bounds `[zMin zMax]`

## Returns
+ `gram`  positive Hilbert-majorant Gram matrix

## Discussion

  The majorant retains the positive interior contribution and
  replaces every signed endpoint coefficient by its absolute
  value. It is the positive product associated with the natural
  $$L^2\oplus\mathbb C^s$$ decomposition of the Pontryagin
  space. It coincides with `gramMatrix` when all endpoint
  coefficients are nonnegative.
