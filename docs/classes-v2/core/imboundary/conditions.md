---
layout: default
title: conditions
parent: IMBoundary
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  conditions

Place bottom and surface boundary conditions for one variable.


---

## Declaration
```matlab
 boundaryConditions = IMBoundary.conditions(options)
```
## Parameters
+ `options.formulation`  EVP formulation, `"F"` or `"G"`
+ `options.surface`  location-free surface boundary law
+ `options.bottom`  location-free bottom boundary law

## Returns
+ `boundaryConditions`  placed boundary-condition array

## Discussion

  The returned array is ordered bottom first and surface second,
  matching the usual vertical domain order.
