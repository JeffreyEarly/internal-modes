---
layout: default
title: GaussQuadraturePointsForEigenmatrices
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 24
mathjax: true
---

#  GaussQuadraturePointsForEigenmatrices

Return quadrature points inferred from a generalized EVP.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 z_g = GaussQuadraturePointsForEigenmatrices(self,nPoints,A,B)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `nPoints`  number of quadrature points requested
+ `A`  left generalized-eigenproblem matrix
+ `B`  right generalized-eigenproblem matrix

## Returns
+ `z_g`  depth locations of the quadrature points

## Discussion
