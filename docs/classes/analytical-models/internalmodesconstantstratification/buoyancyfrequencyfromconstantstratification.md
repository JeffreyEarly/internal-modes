---
layout: default
title: BuoyancyFrequencyFromConstantStratification
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  BuoyancyFrequencyFromConstantStratification

Estimate `N0` and `rho0` from a constant-stratification profile.


---

## Declaration
```matlab
 [N0,rho0] = BuoyancyFrequencyFromConstantStratification(rho,z_in)
```
## Parameters
+ `rho`  density profile as gridded values, a spline, or a function handle
+ `z_in`  depth grid or domain bounds associated with `rho`

## Returns
+ `N0`  estimated constant buoyancy frequency
+ `rho0`  estimated reference surface density

## Discussion


