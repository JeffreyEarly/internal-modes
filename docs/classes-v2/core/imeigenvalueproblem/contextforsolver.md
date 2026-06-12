---
layout: default
title: contextForSolver
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 6
mathjax: true
---

#  contextForSolver

Return the coefficient context for this EVP and solver.


---

## Declaration
```matlab
 context = contextForSolver(evp,solver)
```
## Parameters
+ `solver`  coordinate-aware internal-mode solver

## Returns
+ `context`  framework coefficient context

## Discussion

  The returned struct starts with the solver context, including
  fields such as `N2`, `dzLogN2`, `zDomain`, and `coordinateKind`.
  The EVP then adds physical constants as `f0` and `g`. Operator
  coefficients, boundary rows, inner-product weights, and
  normalization rules read this context but do not own the solver
  discretization.
