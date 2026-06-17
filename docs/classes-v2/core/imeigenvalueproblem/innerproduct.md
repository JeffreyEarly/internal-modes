---
layout: default
title: innerProduct
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 8
mathjax: true
---

#  innerProduct

Return the scalar inner-product recipe.


---

## Declaration
```matlab
 spec = innerProduct(evp)
```
## Returns
+ `spec`  struct with interior and endpoint metric terms

## Discussion

  The canonical solved-variable inner product is
  $$\langle u_i,u_j\rangle_\mu=
  \int_{z_b}^{z_s}r(z)u_i(z)u_j(z)\,dz+
  \sum_{\ell\in S}D_\ell^{-1}L_\ell[u_i]L_\ell[u_j],$$
  where $$S$$ is the set of active, nondegenerate endpoints,
  $$D_\ell=\sigma_\ell(a_\ell d_\ell-b_\ell c_\ell),$$ and
  $$L_\ell[u]=c_\ell u(z_\ell)-d_\ell p(z_\ell)\frac{\partial u}{\partial z}(z_\ell).$$

  The returned struct has fields `variable`, `interiorWeight`,
  `surfaceWeights`, and `bottomWeights`. `interiorWeight`
  stores $$r(z)$$. The endpoint arrays are the same endpoint
  metric terms returned by `endpointWeights("surface")` and
  `endpointWeights("bottom")`; each endpoint weight stores
  `coefficient`, equal to $$D_\ell^{-1}$$, plus `c`, `d`,
  and `location`.
