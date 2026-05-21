---
layout: default
title: waveModesAtWavenumber
parent: IMEigenvalueProblem
grand_parent: Classes
nav_order: 24
mathjax: true
---

#  waveModesAtWavenumber

Create the wave-mode `G` EVP at fixed horizontal wavenumber.


---

## Declaration
```matlab
 evp = IMEigenvalueProblem.waveModesAtWavenumber(options)
```
## Parameters
+ `options.k`  horizontal wavenumber
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  location-free surface boundary law
+ `options.bottomBoundary`  location-free bottom boundary law

## Returns
+ `evp`  fixed-wavenumber wave-mode `G` EVP

## Discussion

  The physical-coordinate strong form is
  $$G_{zz}-K^2G=\lambda(f_0^2-N^2)G/g$$.
