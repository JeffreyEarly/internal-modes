---
layout: default
title: generalizedEnergyMatrix
parent: IMGeostrophicZeroAPVModesBasis
grand_parent: Core
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

Coefficients for endpoints absent from `endpoints` have no
effect because their form matrices are exactly zero.
