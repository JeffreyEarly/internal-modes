---
layout: default
title: gramMatrix
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 4
mathjax: true
---

#  gramMatrix

Sampled modal Gram matrix $$\Gamma$$.

> Developer documentation: this item describes internal implementation details.


---

## Discussion

The transform computes

$$
\Gamma=A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}
=\Phi^\mathsf{T}W\Phi.
$$

It is an $$n_m\times n_m$$ matrix of sampled inner products among
the retained modes. It enters the full definition
$$A_{\mathrm f}=\Gamma^{-1}\Phi^\mathsf{T}W$$ and is compared with
`targetGramMatrix` by `relativeGramOperatorError`.

```matlab
gramDifference = transform.gramMatrix-transform.targetGramMatrix;
```
