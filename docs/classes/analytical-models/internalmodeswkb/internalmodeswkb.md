---
layout: default
title: InternalModesWKB
parent: InternalModesWKB
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  InternalModesWKB

Initialize the WKB approximation solver.


---

## Declaration
```matlab
 im = InternalModesWKB(options)
```
## Parameters
+ `options.rho`  density profile as gridded values, a spline, or a function handle
+ `options.N2`  buoyancy-frequency function handle used instead of `rho`
+ `options.zIn`  input depth grid or domain bounds
+ `options.zOut`  output depth grid
+ `options.latitude`  latitude in degrees
+ `options.rho0`  reference surface density
+ `options.nModes`  optional cap on the number of modes returned
+ `options.nEVP`  spectral initialization resolution
+ `options.rotationRate`  planetary rotation rate in radians per second
+ `options.g`  gravitational acceleration

## Returns
+ `im`  WKB solver instance

## Discussion


