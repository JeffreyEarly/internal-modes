---
layout: default
title: rotateBoundaryDepth
parent: IMGeostrophicZeroAPVModesBasis
grand_parent: Core
nav_order: 17
mathjax: true
---

#  rotateBoundaryDepth

Diagonalize generalized energy relative to endpoint response.


---

## Declaration
```matlab
 basisSet = rotateBoundaryDepth(boundaryModes,options)
```
## Parameters
+ `options.g0`  finite signed surface coefficient
+ `options.gd`  finite signed bottom coefficient

## Returns
+ `basisSet`  boundary-depth rotated basis

## Discussion

  For each wavenumber page, solve

  $$
  \mathsf H_g\mathbf c^a=\gamma_a\mathsf R_B\mathbf c^a,
  \qquad
  (\mathbf c^a)^T\mathsf R_B\mathbf c^b=\delta^{ab},
  \qquad
  h_0^a=2k^2\gamma_a.
  $$

  ```matlab
  depthModes = boundaryModes.rotateBoundaryDepth(g0=-0.035,gd=0.01);
  ```
