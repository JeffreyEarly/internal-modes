---
layout: default
title: normalization
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 10
mathjax: true
---

#  normalization

Name of the normalization captured by this transform.


---

## Discussion

The columns of `inverseMatrix` were sampled using this basis-set
normalization, so modal coefficients are defined relative to the
same scaling. The value is a snapshot taken when the transform was
built; subsequently changing the source basis normalization does
not modify an existing transform.
