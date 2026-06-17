---
layout: default
title: gramMatrix
parent: IMBasisSet
grand_parent: Core
nav_order: 6
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
  $$M_{ij}=\int_{z_a}^{z_b} r(z)u_i(z)u_j(z)\,dz+
  \text{included endpoint terms}.$$
