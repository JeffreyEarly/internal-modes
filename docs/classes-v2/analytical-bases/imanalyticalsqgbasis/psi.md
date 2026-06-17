---
layout: default
title: psi
parent: IMAnalyticalSQGBasis
grand_parent: Analytical bases
nav_order: 6
mathjax: true
---

#  psi

Evaluate exact SQG streamfunction modes.


---

## Declaration
```matlab
 values = psi(sqg,z)
```
## Parameters
+ `z`  physical coordinate

## Returns
+ `values`  SQG streamfunction values

## Discussion

  The returned array has one row per `z` value and one column per
  retained wavenumber in `k`.
