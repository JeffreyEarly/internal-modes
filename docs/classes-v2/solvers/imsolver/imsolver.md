---
layout: default
title: IMSolver
parent: IMSolver
grand_parent: Solvers
nav_order: 1
mathjax: true
---

#  IMSolver

Define the shared protocol for internal-mode solvers.


---

## Declaration
```matlab
 classdef (Abstract) IMSolver
```
## Discussion

  `IMSolver` owns the solver-independent generalized EVP
  workflow. Concrete subclasses provide the native grid, physical
  derivative matrices, boundary rows, and native-mode evaluation.
  Solvers own the numerical medium and discretization. EVPs own the
  physical constants and combine them with solver context during
  assembly.
