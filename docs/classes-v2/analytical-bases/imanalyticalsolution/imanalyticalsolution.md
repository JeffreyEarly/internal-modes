---
layout: default
title: IMAnalyticalSolution
parent: IMAnalyticalSolution
grand_parent: Analytical bases
nav_order: 1
mathjax: true
---

#  IMAnalyticalSolution

Initialize common analytical-stratification state.


---

## Declaration
```matlab
 solution = IMAnalyticalSolution(options)
```
## Parameters
+ `options.zDomain`  physical vertical domain
+ `options.rotationRate`  nonnegative planetary rotation rate in radians per second
+ `options.latitude`  latitude in degrees north, between -90 and 90
+ `options.g`  gravitational acceleration

## Returns
+ `solution`  analytical solution family

## Discussion
