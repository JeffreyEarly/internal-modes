---
layout: default
title: bottomModesAtWavenumber
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 11
mathjax: true
---

#  bottomModesAtWavenumber

Return the analytical bottom SQG mode for constant stratification.


---

## Declaration
```matlab
 psi = bottomModesAtWavenumber(self,k)
```
## Parameters
+ `self`  InternalModesConstantStratification instance
+ `k`  horizontal wavenumber array

## Returns
+ `psi`  bottom SQG mode evaluated on `zOut`

## Discussion

            size(psi) = [size(k); length(z)]
