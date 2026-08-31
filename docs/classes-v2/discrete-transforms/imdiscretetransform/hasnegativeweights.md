---
layout: default
title: hasNegativeWeights
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 5
mathjax: true
---

#  hasNegativeWeights

Whether at least one quadrature weight is negative.


---

## Discussion

  This is `true` when any $$w_i<0$$. A negative weight means the
  quadrature rule is algebraic rather than a positive weighted sum,
  but it does not by itself invalidate the transform. Inspect
  `targetGramIsPositiveDefinite` and the Gram diagnostics to assess
  the resulting metric.
