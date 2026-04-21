---
layout: default
title: InternalModesSpectral
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 38
mathjax: true
---

#  InternalModesSpectral

Initialize the spectral solver on depth coordinates.


---

## Declaration
```matlab
 im = InternalModesSpectral(options)
```
## Parameters
+ `options.rho`  density profile as gridded values, a spline, or a function handle
+ `options.N2`  buoyancy-frequency function handle used instead of `rho`
+ `options.zIn`  input depth grid or domain bounds
+ `options.zOut`  output depth grid
+ `options.latitude`  latitude in degrees
+ `options.rho0`  reference surface density
+ `options.nModes`  optional cap on the number of modes returned
+ `options.nEVP`  number of collocation points in the spectral EVP
+ `options.rotationRate`  planetary rotation rate in radians per second
+ `options.g`  gravitational acceleration

## Returns
+ `im`  spectral solver instance

## Discussion


