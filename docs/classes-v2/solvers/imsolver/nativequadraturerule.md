---
layout: default
title: nativeQuadratureRule
parent: IMSolver
grand_parent: Solvers
nav_order: 13
mathjax: true
---

#  nativeQuadratureRule

Return the solver's native physical quadrature rule.


---

## Declaration
```matlab
 [z,weights] = nativeQuadratureRule(solver,zBounds)
```
## Parameters
+ `zBounds`  physical integration bounds

## Returns
+ `z`  increasing physical quadrature points
+ `weights`  physical-coordinate quadrature weights

## Discussion

The returned points are increasing in physical coordinate and
`weights` are reordered with them. For a WKB-configured
`IMSolverSpectral`, the points are Chebyshev--Lobatto points in
$$x(z)=\int N(z)\,dz$$ and the weights act directly on values
sampled in physical $$z$$.
