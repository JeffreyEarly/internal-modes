---
layout: default
title: nonnegativeConstraint
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 11
mathjax: true
---

#  nonnegativeConstraint

Whether the optimization required nonnegative weights.


---

## Discussion

  When true, the least-squares problem imposed $$w_k\geq0$$ at every
  fixed point. When false, `transform.hasNegativeWeights` reports
  whether the fitted algebraic rule actually contains negative
  weights.
