---
layout: default
title: ProfileWithDoubleExponentialPlusPycnocline
parent: InternalModes
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  ProfileWithDoubleExponentialPlusPycnocline

Build a two-exponential pycnocline profile used by the built-in benchmarks.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [rhoFunc,N2Func,zIn] = ProfileWithDoubleExponentialPlusPycnocline(rho_0,D,delta_p,z_p,z_s,L_s,L_d,N0,Np,Nd)
```
## Parameters
+ `rho_0`  reference surface density
+ `D`  water-column depth
+ `delta_p`  pycnocline thickness scale
+ `z_p`  pycnocline center depth
+ `z_s`  shallow transition depth
+ `L_s`  shallow exponential scale
+ `L_d`  deep exponential scale
+ `N0`  surface buoyancy frequency
+ `Np`  pycnocline buoyancy-frequency scale
+ `Nd`  deep buoyancy-frequency scale

## Returns
+ `rhoFunc`  density function handle
+ `N2Func`  buoyancy-frequency function handle
+ `zIn`  depth domain for the constructed profile

## Discussion


