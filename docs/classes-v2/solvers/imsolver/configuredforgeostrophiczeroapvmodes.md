---
layout: default
title: configuredForGeostrophicZeroAPVModes
parent: IMSolver
grand_parent: Solvers
nav_order: 5
mathjax: true
---

#  configuredForGeostrophicZeroAPVModes

Return a solver configured for geostrophic zero-APV modes.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 solver = configuredForGeostrophicZeroAPVModes(solver,problem)
```
## Parameters
+ `problem`  geostrophic zero-APV problem

## Returns
+ `solver`  configured solver

## Discussion

  Concrete solvers prepare their native grid, coordinate mapping,
  and derivative matrices for the supplied zero-APV problem.
