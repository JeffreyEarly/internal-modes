---
layout: default
title: IMSolverSpectral
parent: IMSolverSpectral
grand_parent: Solvers
nav_order: 1
mathjax: true
---

#  IMSolverSpectral

Create a coordinate-aware spectral solver.


---

## Declaration
```matlab
 solver = IMSolverSpectral(options)
```
## Parameters
+ `options.nEVP`  number of EVP coefficients
+ `options.coordinateKind`  native coordinate kind

## Returns
+ `solver`  initialized spectral solver

## Discussion

The `"z"` coordinate works for any canonical EVP. The `"wkb"`
and `"density"` coordinates use the problem stratification
`N2`, so they can only be configured with `IMInternalModes`
EVPs or `IMGeostrophicZeroAPVModes` problems.
