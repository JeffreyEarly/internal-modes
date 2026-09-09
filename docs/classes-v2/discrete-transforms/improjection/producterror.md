---
layout: default
title: productError
parent: IMProjection
grand_parent: Discrete transforms
nav_order: 16
mathjax: true
---

#  productError

Measure product projection errors with a positive coefficient norm.


---

## Declaration
```matlab
 errors = productError(self,sampledProducts,referenceCoefficients,productNormSquared)
```
## Parameters
+ `sampledProducts`  finite real or complex sample-by-product values
+ `referenceCoefficients`  finite real or complex column-by-product coefficients
+ `productNormSquared`  nonnegative continuous squared product norms

## Returns
+ `errors`  row of relative coefficient-error magnitudes

## Discussion

  Reference coefficients must use the same signed target pairing, column
  coordinates, and normalization. A zero product norm returns zero only when
  sampled values and reference coefficients are both zero, and infinity
  otherwise. This method does not assert reference convergence or coverage.
