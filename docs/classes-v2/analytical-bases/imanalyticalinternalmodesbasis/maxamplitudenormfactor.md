---
layout: default
title: maxAmplitudeNormFactor
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 15
mathjax: true
---

#  maxAmplitudeNormFactor

Return the maximum amplitude of `F` or `G`.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 factor = maxAmplitudeNormFactor(basisSet,iMode,options)
```
## Parameters
+ `iMode`  retained mode index
+ `options.variable`  `"F"` or `"G"`

## Returns
+ `factor`  maximum raw variable amplitude

## Discussion

This developer utility returns the raw maximum-amplitude
scale
$$s_j=\max_z |V_j^{\mathrm{raw}}(z)|$$
for exact analytical modes on the analytical integration
grid.
