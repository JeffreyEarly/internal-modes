---
layout: default
title: leftOperator
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 15
mathjax: true
---

#  leftOperator

Left differential operator.


---

## Discussion

  `leftOperator` contributes the matrix `A` in
  $$Aq=\lambda Bq.$$
  Coefficients are evaluated with the context returned by
  `contextForSolver`, so coefficient functions may use solver-owned
  fields such as `ctx.N2` and EVP-owned fields such as `ctx.g`.
