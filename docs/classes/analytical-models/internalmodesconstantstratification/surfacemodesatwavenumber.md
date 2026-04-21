---
layout: default
title: surfaceModesAtWavenumber
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  surfaceModesAtWavenumber

Return the analytical surface SQG mode for constant stratification.


---

## Declaration
```matlab
 psi = surfaceModesAtWavenumber(self,k)
```
## Parameters
+ `self`  InternalModesConstantStratification instance
+ `k`  horizontal wavenumber array

## Returns
+ `psi`  surface SQG mode evaluated on `zOut`

## Discussion

            size(psi) = [size(k); length(z)]
