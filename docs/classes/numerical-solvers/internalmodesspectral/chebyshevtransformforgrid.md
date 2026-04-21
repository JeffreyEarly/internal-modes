---
layout: default
title: ChebyshevTransformForGrid
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  ChebyshevTransformForGrid

Build a spectral interpolation map from a Lobatto grid to an output grid.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [T,doesOutputGridSpanDomain] = ChebyshevTransformForGrid(lobatto_grid,output_grid)
```
## Parameters
+ `lobatto_grid`  source Chebyshev-Lobatto grid
+ `output_grid`  target grid inside the source bounds

## Returns
+ `T`  transform from Chebyshev coefficients to values on `output_grid`
+ `doesOutputGridSpanDomain`  true when the output grid spans the full Lobatto domain

## Discussion


