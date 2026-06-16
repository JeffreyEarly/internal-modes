---
layout: default
title: hydrostaticGModes
parent: IMInternalModes
grand_parent: Core
nav_order: 11
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
+ `options.surfaceBoundary`  surface boundary condition
+ `options.bottomBoundary`  bottom boundary condition

## Returns
+ `evp`  hydrostatic `G` EVP

## Discussion

  The canonical scalar form is
  $$-\frac{\partial^2 G}{\partial z^2}(z)
  =\lambda\frac{N^2(z)}{g}G(z).$$
  Solved hydrostatic basis sets install the `geostrophic`
  normalization rule and use it by default because they set
  `modeFamily` to `"geostrophic"`. This factory sets
  `parameters.formulation`, `parameters.f0`, and `parameters.g`.

  ```matlab
  evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
  solver = IMSolverSpectral(nEVP=128);
  basisSet = solver.solveEVP(evp,nModes=4);
  G = basisSet.G(z);
  ```
