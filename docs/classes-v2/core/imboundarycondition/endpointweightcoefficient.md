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
  -\left[
  a_i u(z_i)
  -b_i p(z_i)\frac{\partial u}{\partial z}(z_i)
  \right]
  =
  \lambda\left[
  c_i u(z_i)
  -d_i p(z_i)\frac{\partial u}{\partial z}(z_i)
  \right].
  $$

  the stored properties `a`, `b`, `c`, and `d` define

  $$
  D_i=(-1)^{i+1}(a_i d_i-b_i c_i),
  $$

  where Yassin's endpoint indexing makes the sign positive at
  the bottom and negative at the surface. This method returns
  $$D_i^{-1}$$, the scalar coefficient multiplying

  $$
  \left(
  c_i u(z_i)
  -d_i p(z_i)\frac{\partial u}{\partial z}(z_i)
  \right)^2
  $$

  in the endpoint part of the norm. `IMEigenvalueProblem`
  copies this scalar into the `coefficient` field of each
  `endpointWeights` struct.

  This method returns only the scalar coefficient
  $$D_i^{-1}$$. The full endpoint norm term is assembled by
  `IMEigenvalueProblem.endpointWeights` from this coefficient
  and the boundary properties `c` and `d`.
