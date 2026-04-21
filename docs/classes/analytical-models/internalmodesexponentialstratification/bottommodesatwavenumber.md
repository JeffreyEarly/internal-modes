---
layout: default
title: BottomModesAtWavenumber
parent: InternalModesExponentialStratification
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  BottomModesAtWavenumber

Return the analytical bottom SQG mode for exponential stratification.


---

## Declaration
```matlab
 psi = BottomModesAtWavenumber(self,k)
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
