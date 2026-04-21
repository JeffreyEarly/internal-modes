---
layout: default
title: ProfileWithDoubleGaussianPlusPycnocline
parent: InternalModes
grand_parent: Classes
nav_order: 13
mathjax: true
---

#  ProfileWithDoubleGaussianPlusPycnocline

Build a double-Gaussian pycnocline profile used by the built-in benchmarks.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [rhoFunc,N2Func,zIn] = ProfileWithDoubleGaussianPlusPycnocline(rho_0,D,delta_p,z_p,L_s,L_d,N0,Nq,Np,Nd)
```
## Parameters
+ `rho_0`  reference surface density
+ `D`  water-column depth
+ `delta_p`  pycnocline thickness scale
+ `z_p`  pycnocline center depth
+ `L_s`  shallow Gaussian scale
+ `L_d`  deep Gaussian scale
+ `N0`  surface buoyancy frequency
+ `Nq`  mixed-layer buoyancy-frequency scale
+ `Np`  pycnocline buoyancy-frequency scale
+ `Nd`  deep buoyancy-frequency scale

## Returns
+ `rhoFunc`  density function handle
+ `N2Func`  buoyancy-frequency function handle
+ `zIn`  depth domain for the constructed profile

## Discussion


