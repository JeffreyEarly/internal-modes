---
layout: default
title: gramConditionNumber
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 3
mathjax: true
---

#  gramConditionNumber

Two-norm condition number of the sampled Gram matrix.


---

## Discussion

This is $$\kappa_2(\Gamma)$$ for
$$\Gamma=A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}$$. Large values
indicate that the Galerkin normal equations used to construct
`forwardMatrix` are sensitive to perturbations. This diagnostic
depends on both the sampled modes and the metric, unlike
`inverseMatrixConditionNumber`.
