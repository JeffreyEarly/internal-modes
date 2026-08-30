---
layout: default
title: rotateSurfaceBuoyancy
parent: IMGeostrophicZeroAPVModesBasis
grand_parent: Core
nav_order: 18
mathjax: true
---

#  rotateSurfaceBuoyancy

Diagonalize surface buoyancy relative to generalized energy.


---

## Declaration
```matlab
 basisSet = rotateSurfaceBuoyancy(boundaryModes,options)
```
## Parameters
+ `options.g0`  finite signed surface coefficient
+ `options.gd`  finite signed bottom coefficient

## Returns
+ `basisSet`  surface-buoyancy rotated basis

## Discussion

For every wavenumber page, solve

$$
g\mathsf B_0\mathbf c^a=\chi^a\mathsf H_g\mathbf c^a.
$$

The surface-buoyancy-carrying direction is ordered first.
This rotation leaves `h0` empty.

```matlab
surfaceModes = boundaryModes.rotateSurfaceBuoyancy(g0=-0.035,gd=0.01);
```
