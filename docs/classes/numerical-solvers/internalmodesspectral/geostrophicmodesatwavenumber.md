---
layout: default
title: GeostrophicModesAtWavenumber
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 31
mathjax: true
---

#  GeostrophicModesAtWavenumber

Return the geostrophic interior `G`-modes at fixed horizontal wavenumber.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [F,G,h] = GeostrophicModesAtWavenumber(self,k)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `k`  horizontal wavenumber

## Returns
+ `F`  horizontal-velocity mode matrix on `zOut`
+ `G`  vertical-velocity mode matrix on `zOut`
+ `h`  equivalent-depth row vector

## Discussion


