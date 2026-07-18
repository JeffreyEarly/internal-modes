---
layout: default
title: surfacePressureNormFactor
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 24
mathjax: true
---

#  surfacePressureNormFactor

Return the raw surface `F` value.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 factor = surfacePressureNormFactor(basisSet,iMode)
```
## Parameters
+ `iMode`  retained mode index

## Returns
+ `factor`  raw surface-pressure scale factor

## Discussion

This developer utility returns the exact raw surface value
$$s_j=F_j^{\mathrm{raw}}(z_\mathrm{surface}).$$
It gives unit surface `F` value when that value is finite and
nonzero. If the raw value is unavailable, zero, or nonfinite,
it returns `1` so normalization remains well defined.
