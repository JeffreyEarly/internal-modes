---
layout: default
title: IMMeanDensityAnomalyModes
parent: IMMeanDensityAnomalyModes
grand_parent: Core
nav_order: 1
mathjax: true
---

#  IMMeanDensityAnomalyModes

Create a generalized-energy mean-density-anomaly EVP.


---

## Declaration
```matlab
 evp = IMMeanDensityAnomalyModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.g`  gravitational acceleration
+ `options.g0`  signed finite surface acceleration, zero, or positive infinity
+ `options.gd`  signed finite bottom acceleration, zero, or positive infinity

## Returns
+ `evp`  mean-density-anomaly EVP descriptor

## Discussion
