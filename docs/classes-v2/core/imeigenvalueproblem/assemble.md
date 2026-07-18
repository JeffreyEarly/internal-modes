---
layout: default
title: assemble
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 2
mathjax: true
---

#  assemble

Build the canonical matrix pair on a solver grid.

> Developer documentation: this item describes internal implementation details.


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

$$
-\frac{\partial}{\partial z}\left(p(z)\frac{\partial u}{\partial z}\right)
+q(z)u(z)=\lambda r(z)u(z).
$$

The surface and bottom rows are replaced by the endpoint
conditions using endpoint values of `p`, producing the matrix
pencil $$A q = \lambda B q$$. This method is mainly for solver
implementations, diagnostics, and external eigensolver
experiments; ordinary workflows call `solver.solveEVP`.
