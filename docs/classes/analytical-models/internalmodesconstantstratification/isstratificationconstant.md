---
layout: default
title: IsStratificationConstant
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  IsStratificationConstant

Test whether a supplied profile is close to constant stratification.


---

## Declaration
```matlab
 flag = IsStratificationConstant(rho,z_in)
```
## Parameters
+ `rho`  density profile as gridded values, a spline, or a function handle
+ `z_in`  depth grid or domain bounds associated with `rho`

## Returns
+ `flag`  logical scalar indicating whether the profile is approximately constant-stratification

## Discussion


