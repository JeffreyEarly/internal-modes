---
layout: default
title: robinEnergyCoefficient
parent: IMBoundaryCondition
grand_parent: Core
nav_order: 13
mathjax: true
---

#  robinEnergyCoefficient

Return the ordinary Robin endpoint quadratic coefficient.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 beta = robinEnergyCoefficient(boundary,location)
```
## Discussion

For inactive boundary conditions, this is
`(-1)^(i+1)*a/b` with Yassin's `z_1` bottom and `z_2`
surface indexing.
