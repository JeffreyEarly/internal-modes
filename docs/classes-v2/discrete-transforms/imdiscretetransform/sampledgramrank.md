---
layout: default
title: sampledGramRank
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 13
mathjax: true
---

#  sampledGramRank

Numerical rank of the sampled modal Gram matrix.


---

## Discussion

  `sampledGramRank` is the number of singular values of
  $$\Gamma=A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}$$ above the
  scale-aware numerical rank tolerance. A value smaller than the
  number of retained modes means the sampled rule cannot distinguish
  the complete modal family. In that case `forwardMatrix` is formed
  with a pseudoinverse and `relativeGramOperatorError` is `Inf` so an
  automatic retained-band policy can reject the deficient prefix.
