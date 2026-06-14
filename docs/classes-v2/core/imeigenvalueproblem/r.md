---
layout: default
title: r
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 20
mathjax: true
---

#  r

Metric coefficient multiplying the eigenvalue side.


---

## Discussion

  `r` defines the interior metric in $$\lambda r u$$ and in the
  default scalar inner product. It may be a scalar, a vector on the
  solver grid, or a function handle with signature
  `values = r(z,ctx)`.
