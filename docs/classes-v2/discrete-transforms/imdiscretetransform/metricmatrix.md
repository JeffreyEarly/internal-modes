---
layout: default
title: metricMatrix
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 9
mathjax: true
---

#  metricMatrix

Sample-space bilinear-form matrix $$W$$.

> Developer documentation: this item describes internal implementation details.


---

## Discussion

The metric defines the sampled bilinear form

$$
\langle x,y\rangle_W=x^\mathsf{T}Wy.
$$

For transforms built by `IMBasisSet`, its structure is

$$
W=\operatorname{diag}\!\left(r(z_i)\Delta z_i\right)
+W_{\mathrm{endpoint}},
$$

where supported value-only endpoint terms are represented in
$$W_{\mathrm{endpoint}}$$. The matrix is $$n_z\times n_z$$,
symmetric, and may be indefinite.

```matlab
metricSymmetryError = norm(transform.metricMatrix-transform.metricMatrix.',2);
```
