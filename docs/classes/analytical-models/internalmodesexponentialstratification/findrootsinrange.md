---
layout: default
title: FindRootsInRange
parent: InternalModesExponentialStratification
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  FindRootsInRange

Find analytical eigenvalue roots over a bounded search interval.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 r = FindRootsInRange(self,nu,s,bounds,expMinusDOverB)
```
## Parameters
+ `self`  InternalModesExponentialStratification instance
+ `nu`  order function for the Bessel solution
+ `s`  argument function for the Bessel solution
+ `bounds`  two-element search interval
+ `expMinusDOverB`  exponential lower-bound factor `e^{-D/b}`

## Returns
+ `r`  sorted root vector

## Discussion

                    nu(x) is a function of x [used in the solution J_\nu(s)]
  s(x) is a function of x [used in the solution J_\nu(s)]
  bounds is the [xmin xmax] of the region to search for roots
  x_nu is where the solution transitions from big_nu to
  small_nu
