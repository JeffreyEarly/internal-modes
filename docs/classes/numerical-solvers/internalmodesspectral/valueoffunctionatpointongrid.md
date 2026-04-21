---
layout: default
title: ValueOfFunctionAtPointOnGrid
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 47
mathjax: true
---

#  ValueOfFunctionAtPointOnGrid

Evaluate a Chebyshev series at arbitrary points in its physical domain.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 value = ValueOfFunctionAtPointOnGrid(x0,x,func_cheb)
```
## Parameters
+ `x0`  evaluation points
+ `x`  Lobatto grid defining the physical domain
+ `func_cheb`  Chebyshev coefficients

## Returns
+ `value`  series evaluated at `x0`

## Discussion

                We have the Chebyshev coefficents of function func_cheb, defined on grid x, return the value at x0;
