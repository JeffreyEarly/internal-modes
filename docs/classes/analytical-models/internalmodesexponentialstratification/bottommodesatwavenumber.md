---
layout: default
title: bottomModesAtWavenumber
parent: InternalModesExponentialStratification
grand_parent: Classes
nav_order: 15
mathjax: true
---

#  bottomModesAtWavenumber

Return the analytical bottom SQG mode for exponential stratification.


---

## Declaration
```matlab
 psi = bottomModesAtWavenumber(self,k)
```
## Parameters
+ `self`  InternalModesExponentialStratification instance
+ `k`  horizontal wavenumber array

## Returns
+ `psi`  bottom SQG mode evaluated on `zOut`

## Discussion

            Not done in LaCasce 2012, but the calculation is almost
  identical.
  size(psi) = [size(k); length(z)]
