---
layout: default
title: rightOperator
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 23
mathjax: true
---

#  rightOperator

Right differential operator.


---

## Discussion

  `rightOperator` contributes the matrix `B` in the generalized EVP.
  Standard factories use it for the weighted side of the strong form,
  for example `@(z,ctx) (ctx.f0*ctx.f0 - ctx.N2(z))/ctx.g` in the
  fixed-wavenumber wave problem.
