---
layout: default
title: GeostrophicFModesAtWavenumber
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 28
mathjax: true
---

#  GeostrophicFModesAtWavenumber

Return the geostrophic interior `F`-modes at fixed horizontal wavenumber.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [F,G,h] = GeostrophicFModesAtWavenumber(self,k)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `k`  horizontal wavenumber

## Returns
+ `F`  horizontal-velocity mode matrix on `zOut`
+ `G`  vertical-velocity mode matrix on `zOut`
+ `h`  equivalent-depth row vector

## Discussion


