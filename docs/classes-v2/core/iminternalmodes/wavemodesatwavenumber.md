---
layout: default
title: waveModesAtWavenumber
parent: IMInternalModes
grand_parent: Core
nav_order: 17
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
+ `options.surfaceBoundary`  surface boundary condition
+ `options.bottomBoundary`  bottom boundary condition

## Returns
+ `evp`  fixed-wavenumber `G` EVP

## Discussion

  The canonical scalar form is
  $$-\frac{\partial^2 G}{\partial z^2}(z)+K^2G(z)
  =\lambda\frac{N^2(z)-f_0^2}{g}G(z).$$
  Solved fixed-wavenumber basis sets install the `kConstant`
  normalization rule and use it by default.
  This factory adds `parameters.k` and sets
  `parameters.formulation`, `parameters.f0`, and `parameters.g`.
