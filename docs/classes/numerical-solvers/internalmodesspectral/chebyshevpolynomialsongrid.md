---
layout: default
title: ChebyshevPolynomialsOnGrid
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  ChebyshevPolynomialsOnGrid

Evaluate Chebyshev polynomials and derivatives on an arbitrary grid.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [varargout] = ChebyshevPolynomialsOnGrid(x,N_polys)
```
## Parameters
+ `x`  grid on which to evaluate the basis
+ `N_polys`  number of Chebyshev polynomials requested

## Returns
+ `varargout`  basis matrices for the polynomials and requested derivatives

## Discussion

              Chebyshev Polynomials on Grid
  Compute the the first N Chebyshev polynomials and their derivatives for
  an arbitrary grid x.

  x_norm = ChebyshevPolynomialsOnGrid( x ) with exactly one argument, x,
  returns the x normalized to its typical [-1,1] values.

  T = ChebyshevPolynomialsOnGrid( x, N_polys ) returns the first N_poly
  Chebyshev polynomials for an arbitrary grid x.

  [T, T_x, T_xx,...] = ChebyshevPolynomialsOnGrid( x, N_polys ) returns the
  first N_poly Chebyshev polynomials and their derivatives for an arbitrary
  grid x.

  The returned matrices T, T_xx, etc are size(T) = [length(x) N_polys],
  i.e., the polynomials are given column-wise.
