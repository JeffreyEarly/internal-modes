---
layout: default
title: endpointGramTerms
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 9
mathjax: true
---

#  endpointGramTerms

Prepare rank-one endpoint terms for exact Gram matrices.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 terms = endpointGramTerms(basisSet,options)
```
## Parameters
+ `options.variable`  `"F"` or `"G"`
+ `options.zBounds`  integration bounds `[zMin zMax]`
+ `options.normalization`  normalization rule name or enum value
+ `options.useNormalized`  whether returned values use the normalization

## Returns
+ `terms`  struct array with `location`, `coefficient`, `values`, and `kind`

## Discussion

`endpointGramTerms` returns endpoint vectors over all retained
analytical modes. Solved-form endpoint weights use
$$L_\ell[V_j]=c_\ell V_j(z_\ell)-d_\ell p(z_\ell)\frac{\partial V_j}{\partial z}(z_\ell),$$
and contribute
$$M \leftarrow M+\gamma_\ell L_\ell L_\ell^\mathsf{T}.$$
Catalog endpoint value terms use $$V_j(z_\ell)$$ and
contribute
$$M \leftarrow M+\alpha_\ell V_\ell V_\ell^\mathsf{T}.$$
Endpoint terms are omitted when `zBounds` does not include
that endpoint.
