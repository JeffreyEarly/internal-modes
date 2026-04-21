---
layout: default
title: ChebyshevDifferentiationMatrix
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  ChebyshevDifferentiationMatrix

Return the matrix that differentiates Chebyshev coefficients.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 D = ChebyshevDifferentiationMatrix(n)
```
## Parameters
+ `n`  number of coefficients

## Returns
+ `D`  coefficient-space differentiation matrix

## Discussion

            Chebyshev Differentiation Matrix
  This matrix differentiates the first `n` Chebyshev polynomials.
