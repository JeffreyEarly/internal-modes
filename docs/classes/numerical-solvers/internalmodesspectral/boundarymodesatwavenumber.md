---
layout: default
title: boundaryModesAtWavenumber
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 49
mathjax: true
---

#  boundaryModesAtWavenumber

Return either the surface or bottom boundary mode at fixed wavenumber.


---

## Declaration
```matlab
 psi = boundaryModesAtWavenumber(self,k,isSurface)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `k`  horizontal wavenumber array
+ `isSurface`  logical flag selecting the surface mode when true and the bottom mode when false

## Returns
+ `psi`  requested boundary mode evaluated on `zOut`

## Discussion

  This helper estimates the boundary-grid resolution needed to
  resolve the smallest retained boundary mode, then solves the
  manuscript SQG-style boundary-value problem spectrally.

              Estimate the grid resolution necessary to resolve the
  smallest mode.
