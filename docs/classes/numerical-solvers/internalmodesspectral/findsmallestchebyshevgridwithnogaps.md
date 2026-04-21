---
layout: default
title: FindSmallestChebyshevGridWithNoGaps
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 21
mathjax: true
---

#  FindSmallestChebyshevGridWithNoGaps

Find the coarsest Lobatto grid that covers a target grid without gaps.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 z_lobatto_grid = FindSmallestChebyshevGridWithNoGaps(z)
```
## Parameters
+ `z`  target grid

## Returns
+ `z_lobatto_grid`  Lobatto grid whose intervals each contain at most one target point

## Discussion

            Want to create a chebyshev grid that never has two or more point between
  its points. If that makes sense.
