---
layout: default
title: IntegrateChebyshevVector
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 34
mathjax: true
---

#  IntegrateChebyshevVector

Integrate a Chebyshev series in coefficient space.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 v_p = IntegrateChebyshevVector(v)
```
## Parameters
+ `v`  Chebyshev coefficients

## Returns
+ `v_p`  coefficients of the antiderivative with zero value at `x=-1`

## Discussion

            Taken from cumsum as part of chebfun.
  chebfun/@chebtech/cumsum.m

  Copyright 2017 by The University of Oxford and The Chebfun Developers.
  See http://www.chebfun.org/ for Chebfun information.
