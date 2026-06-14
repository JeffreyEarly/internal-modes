---
layout: default
title: contextForSolver
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 4
mathjax: true
---

#  contextForSolver

Return the coefficient context for this EVP and solver.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 context = contextForSolver(evp,solver)
```
## Parameters
+ `solver`  canonical solver

## Returns
+ `context`  coefficient context

## Discussion

  The context begins with `solver.context()`, adds `zDomain`,
  then copies each field of `metadata`. Coefficient handles such
  as `p(z,ctx)`, `q(z,ctx)`, and `r(z,ctx)` receive this struct.
