---
layout: default
title: relativeGramError
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 12
mathjax: true
---

#  relativeGramError

Measure the sampled Gram matrix against its continuous target.


---

## Discussion

Let $$\Gamma$$ be `gramMatrix`, let $$\Gamma_0$$ be
`targetGramMatrix`, and define

$$
S=\operatorname{diag}\!\left(
\left|\operatorname{diag}\Gamma_0\right|^{-1/2}
\right).
$$

The reported error is

$$
\left\|S(\Gamma-\Gamma_0)S\right\|_2.
$$

Zero means the sampled metric reproduces the target modal Gram
matrix exactly. When `targetGramIsPositiveDefinite` is true, this
is the worst-case relative quadratic-form, or Parseval, error over
the retained modal space. For a signed target it remains a useful
magnitude-scaled discrepancy, but not a positive-norm error.
