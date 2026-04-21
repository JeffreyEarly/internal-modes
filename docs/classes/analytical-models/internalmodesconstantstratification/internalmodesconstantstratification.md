---
layout: default
title: InternalModesConstantStratification
parent: InternalModesConstantStratification
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  InternalModesConstantStratification

Initialize the constant-stratification analytical solver.


---

## Declaration
```matlab
 im = InternalModesConstantStratification(options)
```
## Parameters
+ `options.N0`  constant buoyancy frequency in radians per second
+ `options.zIn`  two-element depth domain `[zBottom zSurface]`
+ `options.zOut`  output depth grid
+ `options.latitude`  latitude in degrees
+ `options.rho0`  reference surface density
+ `options.nModes`  optional cap on the number of modes returned
+ `options.rotationRate`  planetary rotation rate in radians per second
+ `options.g`  gravitational acceleration

## Returns
+ `im`  constant-stratification solver instance

## Discussion


