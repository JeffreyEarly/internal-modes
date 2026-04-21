---
layout: default
title: BarotropicMode
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  BarotropicMode

Evaluate a chosen analytical barotropic branch.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [F0,G0,F2,G2,N2G2,uMaxRatio,wMaxRatio,kConstantRatio,omegaConstantRatio] = BarotropicMode(self,solutionType,k_z,h0)
```
## Parameters
+ `self`  InternalModesConstantStratification instance
+ `solutionType`  one of `linear`, `hyperbolic`, or `trig`
+ `k_z`  barotropic vertical wavenumber
+ `h0`  barotropic equivalent depth

## Returns
+ `F0`  barotropic horizontal mode
+ `G0`  barotropic vertical mode
+ `F2`  depth integral of `F0.^2`
+ `G2`  depth integral of `G0.^2`
+ `N2G2`  depth integral of `N2 .* G0.^2`
+ `uMaxRatio`  ratio from the active normalization to `uMax`
+ `wMaxRatio`  ratio from the active normalization to `wMax`
+ `kConstantRatio`  ratio from the active normalization to `kConstant`
+ `omegaConstantRatio`  ratio from the active normalization to `omegaConstant`

## Discussion


