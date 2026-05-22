---
layout: default
title: hydrostaticFModes
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 12
mathjax: true
---

#  hydrostaticFModes

Create the geostrophic hydrostatic `F`-mode EVP.


---

## Declaration
```matlab
 evp = IMEigenvalueProblem.hydrostaticFModes(options)
```
## Parameters
+ `options.g`  gravitational acceleration in meters per second squared
+ `options.surfaceBoundary`  location-free surface boundary law
+ `options.bottomBoundary`  location-free bottom boundary law

## Returns
+ `evp`  hydrostatic `F` EVP

## Discussion

  The physical-coordinate strong form is
  $$F_{zz}-(\partial_z\log N^2)F_z=-\lambda N^2F/g,\qquad h=1/\lambda.$$
  The solved variable is `F`; the linked diagnostic variable is
  $$G=-gN^{-2}F_z.$$
  This EVP declares the barotropic mode,
  $$F_0(z)=1,\qquad G_0(z)=0,$$
  so the barotropic mode is retained before the positive
  baroclinic modes. The default normalization is
  `Normalization.geostrophic`.

  ```matlab
  solver = IMSolverSpectral(N2=@(z) 1e-5*ones(size(z)), zDomain=[-1000 0], nEVP=64);
  evp = IMEigenvalueProblem.hydrostaticFModes();
  basisSet = solver.solveEVP(evp, nModes=4);
  F = basisSet.F(linspace(-1000, 0, 128).');
  ```
