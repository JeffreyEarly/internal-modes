---
layout: default
title: majorantNorm
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 17
mathjax: true
---

#  majorantNorm

Return the positive Hilbert-majorant norm of coefficients.


---

## Declaration
```matlab
 value = majorantNorm(basisSet,coefficients,options)
```
## Parameters
+ `coefficients`  one coefficient per retained mode
+ `options.variable`  `"F"` or `"G"`
+ `options.zBounds`  integration bounds

## Returns
+ `value`  positive scalar norm

## Discussion

  For coefficient vector $$c$$ this returns
  $$\sqrt{c^*M_+c}$$. The quantity
  $$\sqrt{|c^*Mc|}$$ formed from the signed Gram matrix is not
  a norm for arbitrary states.
