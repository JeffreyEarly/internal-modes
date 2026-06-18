---
layout: default
title: IMHydrostaticBoundaryCondition
parent: IMHydrostaticBoundaryCondition
grand_parent: Core
nav_order: 1
mathjax: true
---

#  IMHydrostaticBoundaryCondition

Create a physical hydrostatic endpoint law.


---

## Declaration
```matlab
 law = IMHydrostaticBoundaryCondition(options)
```
## Parameters
+ `options.a`  coefficient multiplying `F` on the \(h^0\) side
+ `options.b`  coefficient multiplying `G` on the \(h^0\) side
+ `options.c`  coefficient multiplying `F` on the \(1/h\) side
+ `options.d`  coefficient multiplying `G` on the \(1/h\) side
+ `options.e`  coefficient multiplying `G` on the \(h\) side

## Returns
+ `law`  hydrostatic endpoint law

## Discussion
