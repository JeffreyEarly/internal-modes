---
layout: default
title: determinant
parent: IMBoundaryCondition
grand_parent: Core
nav_order: 6
mathjax: true
---

#  determinant

Return the signed endpoint determinant.


---

## Declaration
```matlab
 value = determinant(boundary,location)
```
## Parameters
+ `location`  `"surface"` or `"bottom"`

## Returns
+ `value`  signed determinant

## Discussion

  Yassin's endpoint indexing uses `z_1` for the bottom endpoint
  and `z_2` for the surface endpoint, so the sign
  `(-1)^(i+1)` is positive at the bottom and negative at the
  surface.
