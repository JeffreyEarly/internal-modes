---
layout: default
title: BarotropicModeNormalization
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  BarotropicModeNormalization

Return the requested analytical normalization factor for the barotropic branch.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 A = BarotropicModeNormalization(self,norm,solutionType,k_z,h0)
```
## Parameters
+ `self`  InternalModesConstantStratification instance
+ `norm`  normalization convention from `Normalization`
+ `solutionType`  one of `linear`, `hyperbolic`, or `trig`
+ `k_z`  barotropic vertical wavenumber
+ `h0`  barotropic equivalent depth

## Returns
+ `A`  normalization factor

## Discussion


