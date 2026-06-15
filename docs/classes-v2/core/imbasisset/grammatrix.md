---
layout: default
title: gramMatrix
parent: IMBasisSet
grand_parent: Core
nav_order: 7
mathjax: true
---

#  gramMatrix

Return the full-domain scalar Gram matrix.


---

## Declaration
```matlab
 gram = gramMatrix(basisSet)
```
## Returns
+ `gram`  scalar Gram matrix

## Discussion

  For normalized scalar modes, entries are
  $$M_{ij}=\int r u_i u_j\,dz+
  \sum_\ell \gamma_\ell L_\ell[u_i]L_\ell[u_j].$$
