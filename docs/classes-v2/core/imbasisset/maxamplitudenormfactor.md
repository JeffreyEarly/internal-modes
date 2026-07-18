---
layout: default
title: maxAmplitudeNormFactor
parent: IMBasisSet
grand_parent: Core
nav_order: 11
mathjax: true
---

#  maxAmplitudeNormFactor

Return the maximum scalar amplitude.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 factor = maxAmplitudeNormFactor(basisSet,iMode)
```
## Parameters
+ `iMode`  retained mode index

## Returns
+ `factor`  maximum raw scalar amplitude

## Discussion

This developer utility returns the raw maximum-amplitude
scale
$$s_j=\max_z |u_j^{\mathrm{raw}}(z)|$$
on the basis-set integration grid. Normalization rules can
use this factor to make the largest scalar mode amplitude
equal to one.
