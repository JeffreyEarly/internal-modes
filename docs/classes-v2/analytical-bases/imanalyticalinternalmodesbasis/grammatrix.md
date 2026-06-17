---
layout: default
title: gramMatrix
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 11
mathjax: true
---

#  gramMatrix

Return a Gram matrix for exact `F` or `G` modes.


---

## Declaration
```matlab
 gram = gramMatrix(basisSet,options)
```
## Parameters
+ `options.variable`  `"F"` or `"G"`
+ `options.zBounds`  integration bounds

## Returns
+ `gram`  Gram matrix

## Discussion

  For variable $$V$$, endpoint terms are included only when
  `zBounds` reaches the corresponding physical endpoint:
  $$M_{ij}=\int_{z_a}^{z_b} w(z)V_i(z)V_j(z)\,dz+
  \sum_\ell \gamma_\ell L_\ell[V_i]L_\ell[V_j]+
  \sum_\ell \alpha_\ell V_i(z_\ell)V_j(z_\ell).$$
  Use `endpointGramTerms` to inspect the prepared endpoint
  vectors that generate the rank-one endpoint updates.
