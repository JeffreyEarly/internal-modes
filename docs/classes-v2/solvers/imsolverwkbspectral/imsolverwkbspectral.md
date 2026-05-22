---
layout: default
title: IMSolverWKBSpectral
parent: IMSolverWKBSpectral
grand_parent: Solvers
nav_order: 1
mathjax: true
---

#  IMSolverWKBSpectral

Create a WKB-coordinate spectral solver.


---

## Declaration
```matlab
 solver = IMSolverWKBSpectral(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.nEVP`  number of EVP coefficients

## Returns
+ `solver`  initialized WKB spectral solver

## Discussion
