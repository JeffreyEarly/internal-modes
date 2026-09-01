---
layout: default
title: gramMatrix
parent: IMBasisSet
grand_parent: Core
nav_order: 10
mathjax: true
---

#  gramMatrix

Return a scalar Gram matrix.


---

## Declaration
```matlab
 gram = gramMatrix(basisSet,options)
```
## Parameters
+ `options.zBounds`  integration bounds `[zMin zMax]`

## Returns
+ `gram`  scalar Gram matrix

## Discussion

  With no arguments this uses the full basis-set domain. Passing
  `zBounds=[zMin zMax]` restricts the interior integral to that
  interval and includes endpoint terms only when the interval
  contains the corresponding physical endpoint. For normalized
  scalar modes,

  $$
  M_{ij}(z_a,z_b)=
  \int_{z_a}^{z_b}r(z)u_i(z)u_j(z)\,dz+
  \sum_{\ell\in S,\ z_\ell\in[z_a,z_b]}D_\ell^{-1}L_\ell[u_i]L_\ell[u_j],
  $$

  where included endpoint terms use

  $$
  L_\ell[u_j]=c_\ell u_j(z_\ell)-d_\ell p(z_\ell)\frac{\partial u_j}{\partial z}(z_\ell).
  $$

  Use `endpointGramTerms` to inspect the prepared endpoint
  vectors that generate the rank-one endpoint updates.
