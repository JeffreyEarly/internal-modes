---
layout: default
title: RenormalizeForGoodConditioning
parent: InternalModes
grand_parent: Classes
nav_order: 15
mathjax: true
---

#  RenormalizeForGoodConditioning

Renormalize a mode matrix by matching column norms to a common scale.


---

## Declaration
```matlab
 [G_tilde,gamma] = RenormalizeForGoodConditioning(G)
```
## Parameters
+ `G`  mode matrix whose columns are ordered by mode number

## Returns
+ `G_tilde`  renormalized mode matrix
+ `gamma`  column rescaling factors

## Discussion


