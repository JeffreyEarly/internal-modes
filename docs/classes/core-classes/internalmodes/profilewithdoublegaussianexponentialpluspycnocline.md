---
layout: default
title: ProfileWithDoubleGaussianExponentialPlusPycnocline
parent: InternalModes
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  ProfileWithDoubleGaussianExponentialPlusPycnocline

Build a mixed Gaussian-exponential pycnocline profile used by the built-in benchmarks.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [rhoFunc,N2Func,zIn] = ProfileWithDoubleGaussianExponentialPlusPycnocline(rho_0,D,delta_p,z_p,L_s,L_d,z_T,L_deep,N0,Nq,Np)
```
## Parameters
+ `rho_0`  reference surface density
+ `D`  water-column depth
+ `delta_p`  pycnocline thickness scale
+ `z_p`  pycnocline center depth
+ `L_s`  shallow Gaussian scale
+ `L_d`  deep Gaussian scale
+ `z_T`  transition depth to the exponential tail
+ `L_deep`  deep exponential scale
+ `N0`  surface buoyancy frequency
+ `Nq`  mixed-layer buoyancy-frequency scale
+ `Np`  pycnocline buoyancy-frequency scale

## Returns
+ `rhoFunc`  density function handle
+ `N2Func`  buoyancy-frequency function handle
+ `zIn`  depth domain for the constructed profile

## Discussion


