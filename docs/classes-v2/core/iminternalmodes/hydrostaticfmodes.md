---
layout: default
title: hydrostaticFModes
parent: IMInternalModes
grand_parent: Core
nav_order: 10
mathjax: true
---

#  hydrostaticFModes

Create the hydrostatic `F` internal-mode EVP.


---

## Declaration
```matlab
 evp = IMInternalModes.hydrostaticFModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface boundary condition
+ `options.bottomBoundary`  bottom boundary condition

## Returns
+ `evp`  hydrostatic `F` EVP

## Discussion

  The canonical scalar form is
  $$-\frac{\partial}{\partial z}\left(N^{-2}(z)
  \frac{\partial F}{\partial z}(z)\right)
  =\lambda\frac{F(z)}{g}.$$
  The barotropic zero mode is inferred from the canonical left
  problem during mode selection.
  This factory sets `parameters.formulation` and `parameters.g`;
  `parameters.f0` is supplied by the internal-mode constructor
  default.
