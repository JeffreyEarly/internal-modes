---
layout: default
title: waveModesAtWavenumber
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 28
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
+ `options.k`  horizontal wavenumber in radians per meter
+ `options.f0`  Coriolis parameter in radians per second
+ `options.g`  gravitational acceleration in meters per second squared
+ `options.surfaceBoundary`  location-free surface boundary law
+ `options.bottomBoundary`  location-free bottom boundary law

## Returns
+ `evp`  fixed-wavenumber wave-mode `G` EVP

## Discussion

  This factory fixes the horizontal wavenumber `K=options.k` and
  solves the physical-coordinate strong form
  $$G_{zz}-K^2G=\lambda(f_0^2-N^2)G/g,\qquad h=1/\lambda.$$
  The solved variable is `G`; the linked diagnostic variable is
  $$F=hG_z.$$
  The factory stores `parameters.k`, uses the default
  `Normalization.kConstant` normalization, and places the supplied
  location-free boundary laws on the surface and bottom. Omitted
  boundaries are rigid `G=0` boundaries.

  ```matlab
  solver = IMSolverSpectral(N2=@(z) 1e-5*ones(size(z)), zDomain=[-500 0], nEVP=48);
  evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=1e-4);
  basisSet = solver.solveEVP(evp, nModes=6);
  h = basisSet.h;
  modeNumber = basisSet.modeNumber;
  ```
