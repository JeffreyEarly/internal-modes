---
layout: default
title: rotateBoundaryDepth
parent: IMAnalyticalGeostrophicZeroAPVModesBasis
grand_parent: Analytical bases
nav_order: 17
mathjax: true
---

#  rotateBoundaryDepth

Diagonalize generalized energy relative to endpoint response.


---

## Declaration
```matlab
 basisSet = rotateBoundaryDepth(exactModes,options)
```
## Parameters
+ `options.g0`  finite signed surface coefficient
+ `options.gd`  finite signed bottom coefficient

## Returns
+ `basisSet`  boundary-depth rotated exact basis

## Discussion

  For each wavenumber page,

  $$
  \mathsf H_g\mathbf c^a=\gamma_a\mathsf R_B\mathbf c^a,
  \qquad h_0^a=2k^2\gamma_a.
  $$
