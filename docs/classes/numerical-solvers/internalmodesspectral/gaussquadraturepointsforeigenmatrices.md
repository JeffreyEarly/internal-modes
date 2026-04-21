---
layout: default
title: GaussQuadraturePointsForEigenmatrices
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 26
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

                  Now we just need to find the roots of the n+1 mode.
  For constant stratification this should give back the
  standard Fourier modes, i.e., an evenly spaced grid.

  Note that if the boundary conditions are such that G(0)=0 and
  G(-D)=0, then those two points do not encode any information.
  As such, only the first (nPoints-2) modes will encode any
  useful information. So we'd expect cond(G(:,1:(nPoints-2))))
  to be good (low), but not the next.
