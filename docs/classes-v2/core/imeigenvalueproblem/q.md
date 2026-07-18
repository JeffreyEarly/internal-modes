---
layout: default
title: q
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 16
mathjax: true
---

#  q

Coefficient multiplying the solved variable on the left side.


---

## Discussion

`q` defines the multiplication term in $$q u$$. It may be a
scalar, a vector on the solver grid, or a function handle with
signature `values = q(z,ctx)`.
