---
layout: default
title: maxAmplitudeNormFactor
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 11
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
+ `factor`  maximum absolute variable amplitude

## Discussion

This is $$s_j=\max_z |V_j^{\mathrm{raw}}(z)|$$ for the
requested variable on the basis-set integration grid. If
`variable` is omitted, the solved formulation is used.
