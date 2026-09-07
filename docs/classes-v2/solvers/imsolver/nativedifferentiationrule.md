---
layout: default
title: nativeDifferentiationRule
parent: IMSolver
grand_parent: Solvers
nav_order: 13
mathjax: true
---

#  nativeDifferentiationRule

Return one ordered grid for physical quadrature and differentiation.


---

## Declaration
```matlab
 [z,weights,Dz] = nativeDifferentiationRule(solver,zBounds,derivativeOrder)
```
## Parameters
+ `zBounds`  physical integration bounds
+ `derivativeOrder`  physical derivative order; default `1`

## Returns
+ `z`  increasing physical grid points
+ `weights`  physical-coordinate quadrature weights
+ `Dz`  matrix mapping values on `z` to physical derivatives on `z`

## Discussion

  The returned points increase in physical coordinate. `weights`
  integrate values in that order, and `Dz*values` differentiates
  columns sampled in the same order. This keeps the grid, weights,
  and differentiation matrix as one inseparable numerical rule.
