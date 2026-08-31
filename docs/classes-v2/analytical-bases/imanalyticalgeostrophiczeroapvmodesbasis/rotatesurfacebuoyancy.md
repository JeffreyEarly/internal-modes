---
layout: default
title: rotateSurfaceBuoyancy
parent: IMAnalyticalGeostrophicZeroAPVModesBasis
grand_parent: Analytical bases
nav_order: 18
mathjax: true
---

#  rotateSurfaceBuoyancy

Diagonalize surface buoyancy relative to generalized energy.


---

## Declaration
```matlab
 basisSet = rotateSurfaceBuoyancy(exactModes,options)
```
## Parameters
+ `options.g0`  finite signed surface coefficient
+ `options.gd`  finite signed bottom coefficient

## Returns
+ `basisSet`  surface-buoyancy rotated exact basis

## Discussion

  The surface-carrying direction is ordered first. This rotation
  reports its eigenvalues through `rotationEigenvalues` and leaves
  `h0` empty.
