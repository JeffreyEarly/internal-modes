---
layout: default
title: InternalModesWKBHydrostatic
parent: InternalModesWKBHydrostatic
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  InternalModesWKBHydrostatic

Initialize the hydrostatic WKB approximation.


---

## Declaration
```matlab
 im = InternalModesWKBHydrostatic(options)
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
+ `im`  hydrostatic WKB solver instance

## Discussion


