---
layout: default
title: hydrostaticGModes
parent: IMEigenvalueProblem
grand_parent: Classes
nav_order: 12
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
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  location-free surface boundary law
+ `options.bottomBoundary`  location-free bottom boundary law

## Returns
+ `evp`  zero-frequency hydrostatic `G` EVP

## Discussion

  Hydrostatic `G` modes satisfy
  $$G_{zz}=-\lambda N^2G/g$$ and have no nontrivial null `G`
  mode.
