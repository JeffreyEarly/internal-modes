---
layout: default
title: IMExponentialStratificationSolution
parent: IMExponentialStratificationSolution
grand_parent: Analytical bases
nav_order: 1
mathjax: true
---

#  IMExponentialStratificationSolution

Create an exponential-stratification analytical solution family.


---

## Declaration
```matlab
 solution = IMExponentialStratificationSolution(options)
```
## Parameters
+ `options.N0`  surface buoyancy frequency
+ `options.b`  exponential e-folding depth
+ `options.zDomain`  physical vertical domain with surface at zero
+ `options.rotationRate`  nonnegative planetary rotation rate in radians per second
+ `options.latitude`  latitude in degrees north, between -90 and 90
+ `options.g`  gravitational acceleration

## Returns
+ `solution`  analytical solution family

## Discussion
