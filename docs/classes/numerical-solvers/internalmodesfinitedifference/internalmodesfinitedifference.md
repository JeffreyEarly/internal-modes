---
layout: default
title: InternalModesFiniteDifference
parent: InternalModesFiniteDifference
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  InternalModesFiniteDifference

Initialize the finite-difference solver.


---

## Declaration
```matlab
 im = InternalModesFiniteDifference(options)
```
## Parameters
+ `options.rho`  density profile as gridded values, a spline, or a function handle
+ `options.N2`  buoyancy-frequency function handle used instead of `rho`
+ `options.zIn`  input depth grid or domain bounds
+ `options.zOut`  output depth grid
+ `options.latitude`  latitude in degrees
+ `options.rho0`  reference surface density
+ `options.nModes`  optional cap on the number of modes returned
+ `options.orderOfAccuracy`  formal order of accuracy for the finite-difference stencils
+ `options.rotationRate`  planetary rotation rate in radians per second
+ `options.g`  gravitational acceleration

## Returns
+ `im`  finite-difference solver instance

## Discussion


