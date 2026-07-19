---
layout: default
title: relativeGramOperatorError
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 11
mathjax: true
---

#  relativeGramOperatorError

Measure the worst normalized Gram distortion.


---

## Discussion

Let $$\Gamma$$ be `gramMatrix`, let $$\Gamma_0$$ be
`targetGramMatrix`, and define

$$
S=\operatorname{diag}\!\left(
\left|\operatorname{diag}\Gamma_0\right|^{-1/2}
\right).
$$

With

$$
E=S(\Gamma-\Gamma_0)S,
$$

the reported error is

$$
\left\|S(\Gamma-\Gamma_0)S\right\|_2.
$$

Zero means the sampled metric reproduces every retained modal
inner product exactly. The operator norm asks for the single
normalized combination of retained modes whose Gram value is most
distorted. This differs from the default quadrature-fit
`residualNorm`, which uses the Frobenius norm $$\|E\|_{\mathrm F}$$
and aggregates errors over all individual mode pairs. When
`targetGramIsPositiveDefinite` is true, this is the worst-case
relative quadratic-form, or Parseval, error. For a signed target
it remains a useful magnitude-scaled discrepancy, but not a
positive-norm error.

```matlab
operatorError = transform.relativeGramOperatorError;
```
