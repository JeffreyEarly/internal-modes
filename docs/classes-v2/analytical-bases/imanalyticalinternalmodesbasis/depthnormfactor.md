---
layout: default
title: depthNormFactor
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 7
mathjax: true
---

#  depthNormFactor

Return the volume-only depth normalization factor.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 factor = depthNormFactor(basisSet,iMode)
```
## Parameters
+ `iMode`  retained mode index

## Returns
+ `factor`  positive volume root-mean-square `F` factor

## Discussion

  This developer utility evaluates the positive factor
  $$s_j=\sqrt{D^{-1}\int_{z_b}^{z_s}(F_j^{\mathrm{raw}})^2\,dz}$$
  without endpoint terms. The same factor scales exact `F` and
  `G` modes.
