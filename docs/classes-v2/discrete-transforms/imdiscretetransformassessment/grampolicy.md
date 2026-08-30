---
layout: default
title: gramPolicy
parent: IMDiscreteTransformAssessment
grand_parent: Discrete transforms
nav_order: 6
mathjax: true
---

#  gramPolicy

Normalized-Gram retained-band policy result.


---

## Discussion

This policy is always enabled. Its shipped default tolerance is
$$10^{-2}$$, selected by constant/exponential coordinate sweeps;
callers can override it through `discreteTransform`.
