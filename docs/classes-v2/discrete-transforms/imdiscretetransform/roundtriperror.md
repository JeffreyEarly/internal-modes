---
layout: default
title: roundTripError
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 12
mathjax: true
---

#  roundTripError

Measure recovery of retained modal coefficients.


---

## Discussion

The reported value is

$$
\left\|A_{\mathrm f}A_{\mathrm i}-I_{n_m}\right\|_2.
$$

A small value means coefficients transformed back to sample space
and then forward are recovered accurately. This is an algebraic
consistency check on the transform matrices, not a quadrature
accuracy metric. It can be near roundoff even when
`relativeGramOperatorError` is appreciable. It also does not
measure the reconstruction error of an arbitrary sampled profile;
that profile is generally projected by
$$A_{\mathrm i}A_{\mathrm f}$$.
