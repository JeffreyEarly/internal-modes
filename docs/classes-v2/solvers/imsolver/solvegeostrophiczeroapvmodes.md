---
layout: default
title: solveGeostrophicZeroAPVModes
parent: IMSolver
grand_parent: Solvers
nav_order: 15
mathjax: true
---

#  solveGeostrophicZeroAPVModes

Solve canonical geostrophic zero-APV boundary modes.


---

## Declaration
```matlab
 basisSet = solveGeostrophicZeroAPVModes(solver,problem)
```
## Parameters
+ `problem`  geostrophic zero-APV problem

## Returns
+ `basisSet`  canonical boundary-normalized basis

## Discussion

The operator is factored once for each distinct requested
horizontal wavenumber. Every requested unit endpoint response is solved
in the same multiple-right-hand-side operation. No
generalized-energy coefficient or rotation enters this solve.
