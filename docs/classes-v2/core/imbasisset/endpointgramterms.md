---
layout: default
title: endpointGramTerms
parent: IMBasisSet
grand_parent: Core
nav_order: 6
mathjax: true
---

#  endpointGramTerms

Prepare rank-one endpoint terms for scalar Gram matrices.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 terms = endpointGramTerms(basisSet,options)
```
## Parameters
+ `options.zBounds`  integration bounds `[zMin zMax]`
+ `options.normalization`  normalization rule name
+ `options.useNormalized`  whether returned values use the normalization

## Returns
+ `terms`  struct array with `location`, `coefficient`, `values`, and `kind`

## Discussion

  `endpointGramTerms` returns the endpoint pieces that enter the
  scalar Gram matrix. Each returned element contains a full
  row vector over retained modes, not one pairwise matrix
  contribution. Canonical endpoint weights use

  $$
  L_\ell[u_j]=c_\ell u_j(z_\ell)-d_\ell p(z_\ell)\frac{\partial u_j}{\partial z}(z_\ell),
  $$

  and the Gram matrix applies the rank-one update

  $$
  M \leftarrow M+D_\ell^{-1}L_\ell L_\ell^\mathsf{T}.
  $$

  Endpoint terms are omitted when `zBounds` does not include
  that endpoint.
