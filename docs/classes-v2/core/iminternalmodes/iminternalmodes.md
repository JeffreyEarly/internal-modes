---
layout: default
title: IMInternalModes
parent: IMInternalModes
grand_parent: Core
nav_order: 3
mathjax: true
---

#  IMInternalModes

Create an internal-mode canonical EVP.


---

## Declaration
```matlab
 evp = IMInternalModes(options)
```
## Parameters
+ `options.name`  short EVP name
+ `options.zDomain`  physical vertical domain
+ `options.N2`  buoyancy frequency squared function
+ `options.formulation`  solved variable, `"F"` or `"G"`
+ `options.p`  canonical derivative-flux coefficient
+ `options.q`  canonical left-side value coefficient
+ `options.r`  canonical metric coefficient
+ `options.surfaceBoundary`  surface endpoint condition
+ `options.bottomBoundary`  bottom endpoint condition
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.hFromEigenvalue`  equivalent-depth conversion
+ `options.FfromGz`  diagnostic relation from `G` derivative to `F`
+ `options.GfromFz`  diagnostic relation from `F` derivative to `G`
+ `options.parameters`  named coefficient parameters

## Returns
+ `evp`  internal-mode EVP descriptor

## Discussion

  The constructor copies `options.parameters`, then writes the
  internal-mode fields `f0`, `g`, and `formulation` into the
  parameter struct. These constructor-owned fields override
  same-named entries in `options.parameters`.
