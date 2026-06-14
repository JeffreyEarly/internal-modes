---
layout: default
title: assemble
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 2
mathjax: true
---

#  assemble

the canonical matrix pair on a solver grid.


---

## Declaration
```matlab
 [A,B] = assemble(evp,solver)
```
## Parameters
+ `solver`  canonical EVP solver

## Returns
+ `A`  left matrix
+ `B`  right matrix

## Discussion

  Interior rows discretize
  $$-(p u')' + q u = \lambda r u.$$
  The surface and bottom rows are replaced by the endpoint
  conditions using endpoint values of `p`.
