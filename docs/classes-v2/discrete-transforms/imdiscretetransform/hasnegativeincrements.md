---
layout: default
title: hasNegativeIncrements
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 5
mathjax: true
---

#  hasNegativeIncrements

Whether at least one quadrature increment is negative.


---

## Discussion

This is `true` when any $$\Delta z_i<0$$. A negative increment is
useful to flag because it means the quadrature rule is algebraic
rather than a positive weighted sum, but it does not by itself
invalidate the transform. Inspect `targetGramIsPositiveDefinite`
and the Gram diagnostics to assess the resulting metric.
