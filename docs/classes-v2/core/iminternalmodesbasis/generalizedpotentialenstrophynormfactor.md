---
layout: default
title: generalizedPotentialEnstrophyNormFactor
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 11
mathjax: true
---

#  generalizedPotentialEnstrophyNormFactor

Return the depth-mean generalized-potential-enstrophy factor.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 factor = generalizedPotentialEnstrophyNormFactor(basisSet,iMode)
```
## Parameters
+ `iMode`  retained mode index

## Returns
+ `factor`  generalized-potential-enstrophy normalization factor

## Discussion

  The solved `F` inner product contains the interior and active
  endpoint terms of the generalized potential enstrophy. This
  factor scales that complete inner product by the physical
  depth, so normalized modes satisfy
  $$D^{-1}\langle F_i,F_j\rangle_\alpha=\delta_{ij}.$$
