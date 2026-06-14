---
layout: default
title: waveModesAtWavenumber
parent: IMInternalModes
grand_parent: Core
nav_order: 11
mathjax: true
---

#  waveModesAtWavenumber

Create the fixed-wavenumber wave-mode EVP.


---

## Declaration
```matlab
 evp = IMInternalModes.waveModesAtWavenumber(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.k`  horizontal wavenumber
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface endpoint condition
+ `options.bottomBoundary`  bottom endpoint condition

## Returns
+ `evp`  fixed-wavenumber `G` EVP

## Discussion

  The canonical scalar form is
  $$-G''+K^2G=\lambda(N^2-f_0^2)G/g.$$
  The default normalization is `Normalization.kConstant`, and
  parameters include `k`, `formulation`, `f0`, and `g`.
