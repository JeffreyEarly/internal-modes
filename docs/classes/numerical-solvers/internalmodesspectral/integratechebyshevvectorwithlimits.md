---
layout: default
title: IntegrateChebyshevVectorWithLimits
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 37
mathjax: true
---

#  IntegrateChebyshevVectorWithLimits

Integrate a Chebyshev series between two physical limits.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 s = IntegrateChebyshevVectorWithLimits(v,x,a,b)
```
## Parameters
+ `v`  Chebyshev coefficients
+ `x`  Lobatto grid defining the physical interval
+ `a`  lower integration limit
+ `b`  upper integration limit

## Returns
+ `s`  definite integral between `a` and `b`

## Discussion

                  v are the coefficients of the chebyshev polynomials
  x is the domain, a Gauss-Lobatto grid
  a and b are the lower and upper limits, respectively
