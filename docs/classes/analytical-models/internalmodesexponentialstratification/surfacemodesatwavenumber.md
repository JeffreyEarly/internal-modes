---
layout: default
title: surfaceModesAtWavenumber
parent: InternalModesExponentialStratification
grand_parent: Classes
nav_order: 20
mathjax: true
---

#  surfaceModesAtWavenumber

Return the analytical surface SQG mode for exponential stratification.


---

## Declaration
```matlab
 psi = surfaceModesAtWavenumber(self,k)
```
## Parameters
+ `self`  InternalModesExponentialStratification instance
+ `k`  horizontal wavenumber array

## Returns
+ `psi`  surface SQG mode evaluated on `zOut`

## Discussion

            See LaCasce 2012.
  size(psi) = [size(k); length(z)]
