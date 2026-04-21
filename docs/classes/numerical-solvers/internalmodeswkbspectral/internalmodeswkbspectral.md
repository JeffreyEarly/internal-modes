---
layout: default
title: InternalModesWKBSpectral
parent: InternalModesWKBSpectral
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  InternalModesWKBSpectral

Initialize the WKB-coordinate spectral solver.


---

## Declaration
```matlab
 im = InternalModesWKBSpectral(options)
```
## Parameters
+ `options.rho`  density profile as gridded values, a spline, or a function handle
+ `options.N2`  buoyancy-frequency function handle used instead of `rho`
+ `options.zIn`  input depth grid or domain bounds
+ `options.zOut`  output depth grid
+ `options.latitude`  latitude in degrees
+ `options.rho0`  reference surface density
+ `options.nModes`  optional cap on the number of modes returned
+ `options.nEVP`  number of collocation points in the stretched-coordinate EVP
+ `options.rotationRate`  planetary rotation rate in radians per second
+ `options.g`  gravitational acceleration

## Returns
+ `im`  WKB-coordinate spectral solver instance

## Discussion

  This constructor intentionally preserves the public
  name-value API used by `wave-vortex-model`.


