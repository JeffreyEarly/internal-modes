---
layout: default
title: boundaryModesAtWavenumber
parent: InternalModesFiniteDifference
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  boundaryModesAtWavenumber

Return either the surface or bottom SQG mode at fixed horizontal wavenumber.


---

## Declaration
```matlab
 psi = boundaryModesAtWavenumber(self,k,isSurface)
```
## Parameters
+ `self`  InternalModesFiniteDifference instance
+ `k`  horizontal wavenumber array
+ `isSurface`  logical flag selecting the surface mode when true and the bottom mode when false

## Returns
+ `psi`  requested boundary mode evaluated on `zOut`

## Discussion


