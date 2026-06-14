---
layout: default
title: IMSolver
parent: IMSolver
grand_parent: Solvers
nav_order: 1
mathjax: true
---

#  IMSolver

Define the shared protocol for canonical EVP solvers.


---

## Declaration
```matlab
 classdef (Abstract) IMSolver
```
## Discussion

  Concrete solvers own the grid, coordinate mapping, derivative
  matrices, integration rule, and interpolation of native modes. The
  base class owns the common generalized-eigenvalue workflow.
