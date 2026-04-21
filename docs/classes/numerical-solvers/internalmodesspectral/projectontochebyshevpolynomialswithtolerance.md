---
layout: default
title: ProjectOntoChebyshevPolynomialsWithTolerance
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 42
mathjax: true
---

#  ProjectOntoChebyshevPolynomialsWithTolerance

Project a profile onto Chebyshev coefficients to a requested tolerance.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [zLobatto,rho_zCheb] = ProjectOntoChebyshevPolynomialsWithTolerance(zIn,rhoFunc,tol)
```
## Parameters
+ `zIn`  source depth grid
+ `rhoFunc`  function to sample on candidate Lobatto grids
+ `tol`  coefficient-chopping tolerance

## Returns
+ `zLobatto`  selected Lobatto grid
+ `rho_zCheb`  Chebyshev coefficients on that grid

## Discussion


