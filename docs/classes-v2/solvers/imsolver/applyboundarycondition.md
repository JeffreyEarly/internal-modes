---
layout: default
title: applyBoundaryCondition
parent: IMSolver
grand_parent: Solvers
nav_order: 3
mathjax: true
---

#  applyBoundaryCondition

Apply a placed boundary condition to a matrix pair.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [A,B] = applyBoundaryCondition(solver,A,B,boundaryCondition,options)
```
## Parameters
+ `A`  left EVP matrix
+ `B`  right EVP matrix
+ `boundaryCondition`  placed boundary condition
+ `options.context`  framework coefficient context

## Returns
+ `A`  boundary-conditioned left matrix
+ `B`  boundary-conditioned right matrix

## Discussion

  Active metadata-only boundary conditions do not replace matrix
  rows. Standard placed conditions replace the solver-native row
  associated with their physical location.
