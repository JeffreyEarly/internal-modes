---
layout: default
title: gramMatrix
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 12
mathjax: true
---

#  gramMatrix

Return the signed Gram matrix for `F` or `G`.


---

## Declaration
```matlab
 gram = gramMatrix(basisSet,options)
```
## Parameters
+ `options.variable`  `"F"` or `"G"`
+ `options.zBounds`  integration bounds `[zMin zMax]`

## Returns
+ `gram`  Gram matrix

## Discussion

  With no arguments this uses the solved formulation over the full
  basis-set domain. Use `variable="F"` or `variable="G"` to choose
  a physical variable, and `zBounds=[zMin zMax]` to restrict the
  interior integral. Endpoint terms are included only when the
  interval contains the corresponding physical endpoint:

  $$
  M_{ij}=\int_{z_a}^{z_b} w(z)V_i(z)V_j(z)\,dz+
  \sum_\ell \gamma_\ell L_\ell[V_i]L_\ell[V_j]+
  \sum_\ell \alpha_\ell V_i(z_\ell)V_j(z_\ell).
  $$

  This is the Pontryagin pairing used for modal projection and
  signed invariants. It can be indefinite. Use
  `majorantGramMatrix` when a positive matrix is required for
  magnitudes, error tolerances, or convergence tests.

  Use `endpointGramTerms` to inspect the prepared endpoint
  vectors that generate the rank-one endpoint updates.
  The requested variable must have a known inner product; if it
  does not, this method throws
  `IMInternalModesBasis:UnavailableInnerProduct` rather than
  returning an incomplete Gram matrix.
