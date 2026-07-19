---
layout: default
title: depthError
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 3
mathjax: true
---

#  depthError

Difference between the fitted weight sum and `depthTarget`.


---

## Discussion

The value is

$$
\sum_k w_k-D.
$$

It should be near zero when `depthConstraint` is true.
