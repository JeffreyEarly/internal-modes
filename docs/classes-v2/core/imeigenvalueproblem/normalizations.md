---
layout: default
title: normalizations
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 18
mathjax: true
---

#  normalizations

Named modal normalization rules.


---

## Discussion

  Each field stores a function handle with signature
  `scale = rule(basisSet,iMode)`. The returned scale divides one raw
  mode column after the solver has assembled, selected, and linked the
  retained modes. Factory-created EVPs populate names such as
  `unity`, `kConstant`, `omegaConstant`, `geostrophic`, `wMax`,
  `uMax`, and `surfacePressure` when those rules are meaningful.
