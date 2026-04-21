---
layout: default
title: GeostrophicRigidLidModes
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 35
mathjax: true
---

#  GeostrophicRigidLidModes

Return rigid-lid geostrophic modes with displaced boundaries.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [F,G,h] = GeostrophicRigidLidModes(self,eta0,etad)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `eta0`  imposed surface displacement coefficient
+ `etad`  imposed bottom displacement coefficient

## Returns
+ `F`  horizontal-velocity mode matrix on `zOut`
+ `G`  vertical-velocity mode matrix on `zOut`
+ `h`  equivalent-depth row vector

## Discussion


