---
layout: default
title: assemble
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 2
mathjax: true
---

#  assemble

Build generalized-EVP matrices on a solver's native basis.


---

## Declaration
```matlab
 [A,B] = assemble(evp,solver)
```
## Parameters
+ `solver`  coordinate-aware internal-mode solver

## Returns
+ `A`  left generalized-EVP matrix
+ `B`  right generalized-EVP matrix

## Discussion

  `assemble` evaluates `leftOperator` and `rightOperator` with the
  merged solver/EVP context and returns the matrices for
  $$Aq=\lambda Bq.$$
  Interior rows come from the operator discretization. Boundary
  rows are then replaced by the resolved endpoint laws through
  the solver, so a rigid `G` boundary imposes the endpoint-value row for
  `G=0` while active or free boundaries can also declare endpoint
  contributions used by normalization and mode indexing.
