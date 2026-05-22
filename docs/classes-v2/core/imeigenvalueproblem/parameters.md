---
layout: default
title: parameters
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 19
mathjax: true
---

#  parameters

Stored factory-specific physical inputs.


---

## Discussion

  `parameters` records physical inputs supplied to a standard factory
  but not otherwise stored as first-class EVP properties, such as
  `parameters.k` or `parameters.omega`. The EVP identity is `name`,
  boundary laws live in `boundaryConditions`, and physical constants
  live in `f0` and `g`. Core EVP assembly does not consume this struct.
