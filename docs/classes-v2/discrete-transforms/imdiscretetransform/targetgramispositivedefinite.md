---
layout: default
title: targetGramIsPositiveDefinite
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 13
mathjax: true
---

#  targetGramIsPositiveDefinite

Whether the target modal Gram matrix defines a positive norm.


---

## Discussion

`targetGramMatrix` is required to be diagonal, so this property is
true exactly when every target modal norm
$$C_j=(\Gamma_0)_{jj}$$ is positive. A false value does not make
the transform invalid: canonical endpoint weights can produce a
signed metric. It changes the interpretation of
`relativeGramError` from a relative norm error to a signed,
magnitude-scaled Gram discrepancy.
