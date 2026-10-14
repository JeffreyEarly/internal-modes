---
layout: default
title: IMConstantStratificationSolution
parent: IMConstantStratificationSolution
grand_parent: Analytical bases
nav_order: 1
mathjax: true
---

#  IMConstantStratificationSolution

Create a constant-stratification analytical solution family.


---

## Declaration
```matlab
 solution = IMConstantStratificationSolution(options)
```
## Parameters
+ `options.N0`  constant buoyancy frequency
+ `options.zDomain`  physical vertical domain
+ `options.rotationRate`  nonnegative planetary rotation rate in radians per second
+ `options.latitude`  latitude in degrees north, between -90 and 90
+ `options.g`  gravitational acceleration

## Returns
+ `solution`  analytical solution family

## Discussion
