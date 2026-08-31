---
layout: default
title: IsChebyshevGrid
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 37
mathjax: true
---

#  IsChebyshevGrid

Test whether a grid is a Chebyshev-Lobatto grid up to tolerance.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 bool = IsChebyshevGrid(z_in)
```
## Parameters
+ `z_in`  candidate grid

## Returns
+ `bool`  true when `z_in` matches a Lobatto grid

## Discussion

            make sure the grid is monotonically decreasing
