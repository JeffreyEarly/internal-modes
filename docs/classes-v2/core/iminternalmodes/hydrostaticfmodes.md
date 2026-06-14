---
layout: default
title: hydrostaticFModes
parent: IMInternalModes
grand_parent: Core
nav_order: 7
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
+ `options.surfaceBoundary`  surface endpoint condition
+ `options.bottomBoundary`  bottom endpoint condition

## Returns
+ `evp`  hydrostatic `F` EVP

## Discussion

  The canonical scalar form is
  $$-\partial_z(N^{-2}F_z)=\lambda F/g.$$
  It declares the barotropic zero mode.
  Metadata includes `formulation` and `g`.
