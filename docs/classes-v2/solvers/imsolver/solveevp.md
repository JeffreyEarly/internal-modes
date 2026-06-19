---
layout: default
title: solveEVP
parent: IMSolver
grand_parent: Solvers
nav_order: 13
mathjax: true
---

#  solveEVP

Solve an EVP and return a basis set.


---

## Declaration
```matlab
 basisSet = solveEVP(solver,evp,options)
```
## Parameters
+ `evp`  canonical EVP descriptor
+ `options.nModes`  number of modes to retain

## Returns
+ `basisSet`  solved basis set

## Discussion

  If the assembled matrices produce no finite real eigenvalues,
  `solveEVP` throws a matrix-level diagnostic before returning.
