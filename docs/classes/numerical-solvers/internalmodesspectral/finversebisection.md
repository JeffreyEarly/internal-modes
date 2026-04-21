---
layout: default
title: fInverseBisection
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 52
mathjax: true
---

#  fInverseBisection

Invert a monotonic function by vectorized bisection.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 y = fInverseBisection(f,x,yMin,yMax,tol)
```
## Parameters
+ `f`  monotonic function handle to invert
+ `x`  target values
+ `yMin`  lower search bound
+ `yMax`  upper search bound
+ `tol`  termination tolerance

## Returns
+ `y`  approximate inverse values satisfying `f(y) = x`

## Discussion

                   FINVERSEBISECTION(F, X)   Compute F^{-1}(X) using Bisection.
  Taken from cumsum as part of chebfun.
  chebfun/inv.m

  Copyright 2017 by The University of Oxford and The Chebfun Developers.
  See http://www.chebfun.org/ for Chebfun information.
