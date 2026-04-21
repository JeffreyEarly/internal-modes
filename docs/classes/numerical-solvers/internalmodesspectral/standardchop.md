---
layout: default
title: standardChop
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 59
mathjax: true
---

#  standardChop

Choose a truncation index for a Chebyshev coefficient sequence.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 cutoff = standardChop(coeffs,tol)
```
## Parameters
+ `coeffs`  coefficient sequence
+ `tol`  chopping tolerance

## Returns
+ `cutoff`  retained coefficient count

## Discussion

              Copyright 2017 by The University of Oxford and The Chebfun Developers.
  See http://www.chebfun.org/ for Chebfun information.

  This is taken, without comments and safe checks, from the
  above developers. They get sole credit. Jared Aurentz and Nick Trefethen, July 2015.
