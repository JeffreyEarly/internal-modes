---
layout: default
title: InternalModesDensitySpectral
parent: InternalModesDensitySpectral
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  InternalModesDensitySpectral

Initialize the density-coordinate spectral solver.


---

## Declaration
```matlab
 im = InternalModesDensitySpectral(options)
```
## Parameters
+ `options.rho`  density profile as gridded values, a spline, or a function handle
+ `options.N2`  buoyancy-frequency function handle used instead of `rho`
+ `options.zIn`  input depth grid or domain bounds
+ `options.zOut`  output depth grid
+ `options.latitude`  latitude in degrees
+ `options.rho0`  reference surface density
+ `options.nModes`  optional cap on the number of modes returned
+ `options.nEVP`  number of collocation points in the density-coordinate EVP
+ `options.rotationRate`  planetary rotation rate in radians per second
+ `options.g`  gravitational acceleration

## Returns
+ `im`  density-coordinate spectral solver instance

## Discussion


