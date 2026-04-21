---
layout: default
title: BaroclinicModesWithEigenvalue
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  BaroclinicModesWithEigenvalue

Evaluate the analytical baroclinic mode shapes for given eigenvalues.


---

## Declaration
```matlab
 [F,G,F2,G2,N2G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = BaroclinicModesWithEigenvalue(self,k_z,h)
```
## Parameters
+ `self`  InternalModesConstantStratification instance
+ `k_z`  vertical wavenumber row vector
+ `h`  equivalent-depth row vector

## Returns
+ `F`  horizontal-velocity mode matrix
+ `G`  vertical-velocity mode matrix
+ `F2`  depth integrals of `F.^2`
+ `G2`  depth integrals of `G.^2`
+ `N2G2`  depth integrals of `N2 .* G.^2`
+ `uMaxRatio`  ratio from the active normalization to `uMax`
+ `wMaxRatio`  ratio from the active normalization to `wMax`
+ `kConstantRatio`  ratio from the active normalization to `kConstant`
+ `omegaConstantRatio`  ratio from the active normalization to `omegaConstant`

## Discussion


