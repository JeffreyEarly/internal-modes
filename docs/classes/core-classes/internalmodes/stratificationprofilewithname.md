---
layout: default
title: StratificationProfileWithName
parent: InternalModes
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  StratificationProfileWithName

Return one of the built-in benchmark stratification profiles.


---

## Declaration
```matlab
 [rhoFunc,N2Func,zIn] = StratificationProfileWithName(stratification)
```
## Parameters
+ `stratification`  profile name such as `constant`, `exponential`, `pycnocline-constant`, or `latmix-site1`

## Returns
+ `rhoFunc`  density function handle
+ `N2Func`  buoyancy-frequency function handle
+ `zIn`  depth domain or depth grid associated with the profile

## Discussion


