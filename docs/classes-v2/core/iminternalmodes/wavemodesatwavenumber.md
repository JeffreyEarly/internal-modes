---
layout: default
title: waveModesAtWavenumber
parent: IMInternalModes
grand_parent: Core
nav_order: 10
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
  Solved fixed-wavenumber basis sets install the `kConstant`
  normalization rule and use it by default.
  This factory adds `parameters.k` and sets
  `parameters.formulation`, `parameters.f0`, and `parameters.g`.
