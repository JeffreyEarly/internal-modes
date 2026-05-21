---
layout: default
title: contextForSolver
parent: IMEigenvalueProblem
grand_parent: Classes
nav_order: 5
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

  The solver supplies medium and discretization fields; the EVP
  supplies physical constants.
