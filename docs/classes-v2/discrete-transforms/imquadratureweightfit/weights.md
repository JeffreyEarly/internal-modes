---
layout: default
title: weights
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 19
mathjax: true
---

#  weights

Fitted algebraic quadrature weights aligned with the fixed points.


---

## Discussion

  These are the weights $$w_k$$ returned as the first output of
  `quadratureWeightsForPoints`. They minimize the selected objective
  subject to the requested constraints and satisfy
  `weights = transform.weights`.
