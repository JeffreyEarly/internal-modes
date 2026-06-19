---
layout: default
title: solveSurfaceGeostrophicModes
parent: IMSolver
grand_parent: Solvers
nav_order: 14
mathjax: true
---

#  solveSurfaceGeostrophicModes

Solve surface-geostrophic boundary modes.


---

## Declaration
```matlab
 basisSet = solveSurfaceGeostrophicModes(solver,problem)
```
## Parameters
+ `problem`  surface-geostrophic boundary-mode problem

## Returns
+ `basisSet`  solved surface-geostrophic basis

## Discussion

  `solveSurfaceGeostrophicModes` solves the boundary-value
  problem stored by `IMSurfaceGeostrophicModes` and returns an
  `IMSurfaceGeostrophicModesBasis`.
