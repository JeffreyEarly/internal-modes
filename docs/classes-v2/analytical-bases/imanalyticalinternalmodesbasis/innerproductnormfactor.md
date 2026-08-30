---
layout: default
title: innerProductNormFactor
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 14
mathjax: true
---

#  innerProductNormFactor

Return the raw inner-product norm factor.

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
+ `factor`  raw inner-product scale factor

## Discussion

This developer utility returns the raw scale factor
$$s_j=\sqrt{|\langle V_j,V_j\rangle|}$$
for exact analytical `F` or `G` modes before the active
normalization is applied. The requested variable must have a
known inner product.
