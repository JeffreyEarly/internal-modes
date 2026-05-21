---
layout: default
title: hydrostaticFModes
parent: IMEigenvalueProblem
grand_parent: Classes
nav_order: 11
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
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  location-free surface boundary law
+ `options.bottomBoundary`  location-free bottom boundary law

## Returns
+ `evp`  hydrostatic `F` EVP

## Discussion

  The physical-coordinate strong form is
  $$F_{zz}-(\partial_z\log N^2)F_z=-\lambda N^2F/g$$.
