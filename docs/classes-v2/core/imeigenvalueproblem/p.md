---
layout: default
title: p
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 19
mathjax: true
---

#  p

Coefficient multiplying the derivative flux.


---

## Discussion

  `p` defines the flux term in $$-(p u')'$$. It may be a scalar,
  a vector on the solver grid, or a function handle with signature
  `values = p(z,ctx)`.
