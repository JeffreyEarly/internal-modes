---
layout: default
title: solveSurfaceGeostrophicModes
parent: IMSolver
grand_parent: Solvers
nav_order: 14
mathjax: true
---

#  solveSurfaceGeostrophicModes

Solve projected surface-geostrophic boundary modes.


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

  `solveSurfaceGeostrophicModes` solves the raw zero-APV
  endpoint modes stored by `IMSurfaceGeostrophicModes`, forms
  the boundary-energy projection, and returns an
  `IMSurfaceGeostrophicModesBasis` with `F`, `G`, and `h`.
