---
layout: default
title: weights
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 17
mathjax: true
---

#  weights

Quadrature weights associated with the sample points.


---

## Discussion

`weights(i)` is the quadrature weight $$w_i$$ associated with
`z(i)`. For transforms built by `IMBasisSet`, the interior part of
the sampled metric begins with

$$
W_{\mathrm{int}}
=\operatorname{diag}\!\left(r(z_i)w_i\right),
$$

before supported endpoint terms are added. The weights may be
geometric control-volume widths or fitted algebraic weights.
`hasNegativeWeights` reports whether any are negative.
