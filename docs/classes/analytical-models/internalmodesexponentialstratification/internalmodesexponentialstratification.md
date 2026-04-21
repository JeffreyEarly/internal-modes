---
layout: default
title: InternalModesExponentialStratification
parent: InternalModesExponentialStratification
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  InternalModesExponentialStratification

Initialize the exponential-stratification analytical solver.


---

## Declaration
```matlab
 im = InternalModesExponentialStratification(options)
```
## Parameters
+ `options.N0`  surface buoyancy frequency in radians per second
+ `options.b`  exponential e-folding scale in meters
+ `options.zIn`  two-element depth domain `[zBottom zSurface]`
+ `options.zOut`  output depth grid
+ `options.latitude`  latitude in degrees
+ `options.rho0`  reference surface density
+ `options.nModes`  optional cap on the number of modes returned
+ `options.rotationRate`  planetary rotation rate in radians per second
+ `options.g`  gravitational acceleration

## Returns
+ `im`  exponential-stratification solver instance

## Discussion


