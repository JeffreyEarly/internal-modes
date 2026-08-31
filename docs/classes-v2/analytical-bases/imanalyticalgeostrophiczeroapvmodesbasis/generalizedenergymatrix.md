---
layout: default
title: generalizedEnergyMatrix
parent: IMAnalyticalGeostrophicZeroAPVModesBasis
grand_parent: Analytical bases
nav_order: 11
mathjax: true
---

#  generalizedEnergyMatrix

Return $$\mathsf H_g=\mathsf H+g_0\mathsf B_0+g_d\mathsf B_d$$.


---

## Declaration
```matlab
 matrix = generalizedEnergyMatrix(basisSet,options)
```
## Parameters
+ `options.g0`  finite signed surface coefficient
+ `options.gd`  finite signed bottom coefficient

## Returns
+ `matrix`  generalized-energy matrix pages

## Discussion

  A coefficient for an absent endpoint has no effect because its
  form matrix is exactly zero.
