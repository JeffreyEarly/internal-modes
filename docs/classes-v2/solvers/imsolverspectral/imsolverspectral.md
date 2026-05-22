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
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.nEVP`  number of EVP coefficients
+ `options.coordinateKind`  native coordinate kind

## Returns
+ `solver`  initialized spectral solver

## Discussion
