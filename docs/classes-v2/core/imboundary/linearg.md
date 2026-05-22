---
layout: default
title: linearG
parent: IMBoundary
grand_parent: Core
nav_order: 21
mathjax: true
---

#  linearG

Create a location-free linear `G` boundary law.


---

## Declaration
```matlab
 boundary = IMBoundary.linearG(options)
```
## Parameters
+ `options.a`  coefficient multiplying `G_z`
+ `options.b`  eigenvalue coefficient multiplying `G`
+ `options.c`  eigenvalue coefficient multiplying `G_z`
+ `options.e`  coefficient multiplying `G`

## Returns
+ `boundary`  initialized linear `G` boundary condition

## Discussion

  When placed, the assembled boundary condition represents
  $$-(eG-aG_z)=\lambda(bG-cG_z)/g$$. Supported pairwise cases
  also declare compatible boundary inner-product terms. General
  unresolved coefficient patterns can still be solved, but
  placing them warns because their inner-product contribution is
  unknown. The derivative $$G_z$$ always means
  $$\partial_z G$$ at either endpoint.
