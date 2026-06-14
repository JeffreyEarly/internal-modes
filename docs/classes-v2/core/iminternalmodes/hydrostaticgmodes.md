---
layout: default
title: hydrostaticGModes
parent: IMInternalModes
grand_parent: Core
nav_order: 9
mathjax: true
---

#  hydrostaticGModes

Create the hydrostatic `G` internal-mode EVP.


---

## Declaration
```matlab
 evp = IMInternalModes.hydrostaticGModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface endpoint condition
+ `options.bottomBoundary`  bottom endpoint condition

## Returns
+ `evp`  hydrostatic `G` EVP

## Discussion

  The canonical scalar form is
  $$-G''=\lambda N^2G/g.$$
  The default normalization is `Normalization.geostrophic`, and
  metadata includes `formulation`, `f0`, and `g`.

  ```matlab
  evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
  solver = IMSolverSpectral(nEVP=128);
  basisSet = solver.solveEVP(evp,nModes=4);
  G = basisSet.G(z);
  ```
