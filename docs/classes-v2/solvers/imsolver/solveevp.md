---
layout: default
title: solveEVP
parent: IMSolver
grand_parent: Solvers
nav_order: 12
mathjax: true
---

#  solveEVP

Solve an EVP and return a native-basis solution set.


---

## Declaration
```matlab
 basisSet = solveEVP(solver,evp,options)
```
## Parameters
+ `evp`  physical-coordinate EVP descriptor
+ `options.nModes`  number of modes to retain

## Returns
+ `basisSet`  solved native basis set

## Discussion

  If the assembled matrices produce no finite real eigenvalues,
  `solveEVP` throws an explanatory diagnostic before returning.
