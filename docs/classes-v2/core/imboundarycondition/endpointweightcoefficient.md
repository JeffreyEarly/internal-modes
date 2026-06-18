---
layout: default
title: endpointWeightCoefficient
parent: IMBoundaryCondition
grand_parent: Core
nav_order: 9
mathjax: true
---

#  endpointWeightCoefficient

Return the endpoint weight coefficient.


---

## Declaration
```matlab
 value = endpointWeightCoefficient(boundary,location)
```
## Parameters
+ `location`  `"surface"` or `"bottom"`

## Returns
+ `value`  coefficient multiplying `(c*u-d*p*u_z)^2`

## Discussion

  For the active boundary condition

  $$
  -[a_i u-b_i(pu')]=\lambda[c_i u-d_i(pu')].
  $$

  the stored properties `a`, `b`, `c`, and `d` define

  $$
  D_i=(-1)^{i+1}(a_i d_i-b_i c_i),
  $$

  where Yassin's endpoint indexing makes the sign positive at
  the bottom and negative at the surface. This method returns
  $$D_i^{-1}$$, the scalar coefficient multiplying

  $$
  (c_i u-d_i p u_z)^2
  $$

  in the endpoint part of the norm. `IMEigenvalueProblem`
  copies this scalar into the `coefficient` field of each
  `endpointWeights` struct.

  This method returns only the scalar coefficient
  $$D_i^{-1}$$. The full endpoint norm term is assembled by
  `IMEigenvalueProblem.endpointWeights` from this coefficient
  and the boundary properties `c` and `d`.
