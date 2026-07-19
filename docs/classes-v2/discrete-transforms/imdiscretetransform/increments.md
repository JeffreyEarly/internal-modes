---
layout: default
title: increments
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 6
mathjax: true
---

#  increments

Quadrature increments associated with the sample points.


---

## Discussion

`increments(i)` is the increment $$\Delta z_i$$ associated with
`z(i)`. For transforms built by `IMBasisSet`, the interior part of
the sampled metric begins with

$$
W_{\mathrm{int}}
=\operatorname{diag}\!\left(r(z_i)\Delta z_i\right),
$$

before supported endpoint terms are added. Signed increments are
retained so algebraic quadrature rules can be inspected directly;
`hasNegativeIncrements` reports whether any are negative.
