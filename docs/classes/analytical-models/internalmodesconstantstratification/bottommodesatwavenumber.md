---
layout: default
title: BottomModesAtWavenumber
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  BottomModesAtWavenumber

Return the analytical bottom SQG mode for constant stratification.


---

## Declaration
```matlab
 psi = BottomModesAtWavenumber(self,k)
```
## Parameters
+ `self`  InternalModesConstantStratification instance
+ `k`  horizontal wavenumber array

## Returns
+ `psi`  bottom SQG mode evaluated on `zOut`

## Discussion

            size(psi) = [size(k); length(z)]
