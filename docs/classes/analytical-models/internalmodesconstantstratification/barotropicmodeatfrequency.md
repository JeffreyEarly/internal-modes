---
layout: default
title: BarotropicModeAtFrequency
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  BarotropicModeAtFrequency

Return the analytical barotropic mode branch for fixed $$\omega$$.


---

## Declaration
```matlab
 [F0,G0,h0,F20,G20,N2G20,uMaxRatio0,wMaxRatio0,kConstantRatio0,omegaConstantRatio0] = BarotropicModeAtFrequency(self,omega)
```
## Parameters
+ `self`  InternalModesConstantStratification instance
+ `omega`  frequency in radians per second

## Returns
+ `F0`  barotropic horizontal mode
+ `G0`  barotropic vertical mode
+ `h0`  barotropic equivalent depth
+ `F20`  depth integral of `F0.^2`
+ `G20`  depth integral of `G0.^2`
+ `N2G20`  depth integral of `N2 .* G0.^2`
+ `uMaxRatio0`  ratio from the active normalization to `uMax`
+ `wMaxRatio0`  ratio from the active normalization to `wMax`
+ `kConstantRatio0`  ratio from the active normalization to `kConstant`
+ `omegaConstantRatio0`  ratio from the active normalization to `omegaConstant`

## Discussion


