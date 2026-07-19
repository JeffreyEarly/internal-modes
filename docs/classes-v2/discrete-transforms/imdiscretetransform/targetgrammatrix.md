---
layout: default
title: targetGramMatrix
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 14
mathjax: true
---

#  targetGramMatrix

Continuous diagonal Gram matrix targeted by the quadrature rule.

> Developer documentation: this item describes internal implementation details.


---

## Discussion

The matrix

$$
\Gamma_0=\operatorname{diag}(C_1,\ldots,C_{n_m})
$$

contains the continuous full-domain inner products of the retained
normalized modes. A fitted quadrature rule seeks to make
`gramMatrix` reproduce this target. The entries $$C_j$$ are finite
and nonzero, but may be negative for a signed canonical metric.

```matlab
targetNorms = diag(transform.targetGramMatrix);
```
