---
layout: default
title: IsStratificationExponential
parent: InternalModesExponentialStratification
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  IsStratificationExponential

Test whether a supplied profile is close to the exponential benchmark.


---

## Declaration
```matlab
 [flag,rho_params] = IsStratificationExponential(rho,z_in)
```
## Parameters
+ `rho`  density profile as gridded values, a spline, or a function handle
+ `z_in`  depth grid or domain bounds associated with `rho`

## Returns
+ `flag`  logical scalar indicating whether the profile matches the exponential benchmark
+ `rho_params`  inferred `[N0 b rho0]` parameter vector

## Discussion


