---
layout: default
title: majorantInnerProduct
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 16
mathjax: true
---

#  majorantInnerProduct

Return the induced positive Hilbert-majorant recipe.


---

## Declaration
```matlab
 spec = majorantInnerProduct(basisSet,options)
```
## Parameters
+ `options.variable`  `"F"` or `"G"`

## Returns
+ `spec`  positive interior and absolute-endpoint recipe

## Discussion

  The recipe retains the positive interior weight and replaces
  every signed endpoint coefficient by its absolute value. Use
  `evp.innerProduct` for the signed Pontryagin recipe used by
  projection and signed invariants.
