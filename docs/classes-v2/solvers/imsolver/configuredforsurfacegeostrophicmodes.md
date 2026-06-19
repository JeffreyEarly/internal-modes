---
layout: default
title: configuredForSurfaceGeostrophicModes
parent: IMSolver
grand_parent: Solvers
nav_order: 5
mathjax: true
---

#  configuredForSurfaceGeostrophicModes

Return a solver configured for surface-geostrophic modes.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 solver = configuredForSurfaceGeostrophicModes(solver,problem)
```
## Parameters
+ `problem`  surface-geostrophic boundary-mode problem

## Returns
+ `solver`  configured solver

## Discussion

  Concrete solvers prepare their native grid, coordinate mapping,
  and derivative matrices for the supplied SQG problem.
