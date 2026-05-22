---
layout: default
title: linearF
parent: IMBoundary
grand_parent: Core
nav_order: 20
mathjax: true
---

#  linearF

Create a location-free linear `F` boundary law.


---

## Declaration
```matlab
 boundary = IMBoundary.linearF(options)
```
## Parameters
+ `options.a`  coefficient multiplying `F`
+ `options.b`  coefficient multiplying `F_z/N2`
+ `options.c`  eigenvalue coefficient multiplying `F`
+ `options.d`  eigenvalue coefficient multiplying `F_z/N2`

## Returns
+ `boundary`  initialized linear `F` boundary condition

## Discussion

  When placed, the assembled boundary condition represents
  $$-(aF-bF_z/N^2)=\lambda(cF-dF_z/N^2)/g$$. Supported pairwise
  cases also declare compatible boundary inner-product terms.
  General unresolved coefficient patterns can still be solved,
  but placing them warns because their inner-product contribution
  is unknown. The derivative $$F_z$$ always means
  $$\partial_z F$$ at either endpoint.
