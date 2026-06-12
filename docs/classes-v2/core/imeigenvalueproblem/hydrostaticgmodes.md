---
layout: default
title: hydrostaticGModes
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 14
mathjax: true
---

#  hydrostaticGModes

Create the hydrostatic `G`-mode EVP.


---

## Declaration
```matlab
 evp = IMEigenvalueProblem.hydrostaticGModes(options)
```
## Parameters
+ `options.f0`  Coriolis parameter in radians per second
+ `options.g`  gravitational acceleration in meters per second squared
+ `options.surfaceBoundary`  location-free surface boundary law
+ `options.bottomBoundary`  location-free bottom boundary law

## Returns
+ `evp`  zero-frequency hydrostatic `G` EVP

## Discussion

  Hydrostatic `G` modes are the zero-frequency wave-mode problem
  written directly as
  $$G_{zz}=-\lambda N^2G/g,\qquad h=1/\lambda.$$
  The solved variable is `G`; the linked diagnostic variable is
  $$F=hG_z.$$
  There is no nontrivial null `G` mode, so retained modes are the
  boundary-index modes declared by the boundary laws followed by
  positive interior baroclinic modes. The default normalization is
  `Normalization.geostrophic`.

  ```matlab
  solver = IMSolverSpectral(N2=@(z) 1e-5*ones(size(z)), zDomain=[-1000 0], nEVP=64);
  evp = IMEigenvalueProblem.hydrostaticGModes();
  basisSet = solver.solveEVP(evp, nModes=4);
  G = basisSet.G(linspace(-1000, 0, 128).');
  ```
