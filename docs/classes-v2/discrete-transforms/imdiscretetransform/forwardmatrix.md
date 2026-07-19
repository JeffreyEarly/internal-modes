---
layout: default
title: forwardMatrix
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 2
mathjax: true
---

#  forwardMatrix

Map sampled profiles to retained modal coefficients.


---

## Discussion

The `forwardMatrix` is the $$n_m\times n_z$$ Galerkin matrix

$$
A_{\mathrm f}
=\left(A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}\right)^{-1}
A_{\mathrm i}^\mathsf{T}W
=\Gamma^{-1}\Phi^\mathsf{T}W,
$$

where $$A_{\mathrm i}=\Phi$$ is `inverseMatrix`, $$W$$ is
`metricMatrix`, and $$\Gamma=\Phi^\mathsf{T}W\Phi$$ is
`gramMatrix`. For sampled profiles $$X$$, the coefficients
$$A=A_{\mathrm f}X$$ satisfy the Galerkin normal equations

$$
\Phi^\mathsf{T}W\left(X-\Phi A\right)=0.
$$

Direct multiplication and `transformForward` are equivalent:

```matlab
coefficientsByMatrix = transform.forwardMatrix*values;
coefficients = transform.transformForward(values);
```
