---
layout: default
title: InternalModesAdaptiveSpectral
parent: InternalModesAdaptiveSpectral
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  InternalModesAdaptiveSpectral

Initialize the adaptive WKB spectral solver.


---

## Declaration
```matlab
 im = InternalModesAdaptiveSpectral(options)
```
## Parameters
+ `options.rho`  density profile as gridded values, a spline, or a function handle
+ `options.N2`  buoyancy-frequency function handle used instead of `rho`
+ `options.zIn`  input depth grid or domain bounds
+ `options.zOut`  output depth grid
+ `options.latitude`  latitude in degrees
+ `options.rho0`  reference surface density
+ `options.nModes`  optional cap on the number of modes returned
+ `options.nEVP`  total number of collocation points across all coupled subproblems
+ `options.rotationRate`  planetary rotation rate in radians per second
+ `options.g`  gravitational acceleration

## Returns
+ `im`  adaptive WKB spectral solver instance

## Discussion


