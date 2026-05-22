---
layout: default
title: hFromEigenvalue
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 10
mathjax: true
---

#  hFromEigenvalue

Equivalent-depth conversion function.


---

## Discussion

  `hFromEigenvalue` is a function handle with signature
  `h = hFromEigenvalue(lambda)`. Standard EVPs use
  $$h_j=1/\lambda_j,$$
  so `lambda` has inverse-depth units in those problems.
