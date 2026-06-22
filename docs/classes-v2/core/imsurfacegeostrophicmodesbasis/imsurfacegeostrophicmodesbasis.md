---
layout: default
title: IMSurfaceGeostrophicModesBasis
parent: IMSurfaceGeostrophicModesBasis
grand_parent: Core
nav_order: 3
mathjax: true
---

#  IMSurfaceGeostrophicModesBasis

Create a solved surface-geostrophic basis.


---

## Declaration
```matlab
 basisSet = IMSurfaceGeostrophicModesBasis(options)
```
## Parameters
+ `options.problem`  surface-geostrophic problem descriptor
+ `options.solver`  configured numerical solver
+ `options.nativeModes`  native projected mode columns
+ `options.k`  mode-aligned horizontal wavenumbers
+ `options.h`  equivalent boundary depths
+ `options.modeNumber`  projected mode labels
+ `options.mixingCoefficients`  raw endpoint-mode coefficients
+ `options.energyEigenvalues`  boundary-energy eigenvalues
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  solved surface-geostrophic basis

## Discussion
