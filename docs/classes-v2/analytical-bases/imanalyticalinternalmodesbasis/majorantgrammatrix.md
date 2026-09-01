---
layout: default
title: majorantGramMatrix
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
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
+ `options.zBounds`  integration bounds

## Returns
+ `gram`  positive Hilbert-majorant Gram matrix

## Discussion

  The majorant retains the positive interior contribution and
  replaces every signed endpoint coefficient by its absolute
  value. It coincides with `gramMatrix` when all endpoint
  coefficients are nonnegative.
