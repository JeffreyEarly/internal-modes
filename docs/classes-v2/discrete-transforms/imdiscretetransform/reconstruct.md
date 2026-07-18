---
layout: default
title: reconstruct
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 13
mathjax: true
---

#  reconstruct

Return sampled profiles reconstructed from retained coefficients.


---

## Declaration
```matlab
 values = reconstruct(transform,coefficients)
```
## Parameters
+ `coefficients`  coefficient arrays with rows aligned to `modeNumber`

## Returns
+ `values`  reconstructed profiles sampled on `z`

## Discussion

This method applies `inverseMatrix`:

$$
\widehat{\mathbf{x}}=A_{\mathrm i}\mathbf{a}.
$$
