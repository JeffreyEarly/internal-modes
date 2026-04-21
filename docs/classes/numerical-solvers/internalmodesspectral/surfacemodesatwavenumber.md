---
layout: default
title: surfaceModesAtWavenumber
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 59
mathjax: true
---

#  surfaceModesAtWavenumber

Return the surface SQG mode at fixed horizontal wavenumber.


---

## Declaration
```matlab
 psi = surfaceModesAtWavenumber(self,k)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `k`  horizontal wavenumber array

## Returns
+ `psi`  surface boundary mode evaluated on `zOut`

## Discussion

  Section 5.6 of Early, Lelong, and Smith (2020) discusses the
  surface-trapped modes computed by this helper.


