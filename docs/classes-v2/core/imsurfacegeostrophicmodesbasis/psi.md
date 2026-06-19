---
layout: default
title: psi
parent: IMSurfaceGeostrophicModesBasis
grand_parent: Core
nav_order: 8
mathjax: true
---

#  psi

Evaluate SQG streamfunction modes.


---

## Declaration
```matlab
 values = psi(basisSet,z)
```
## Parameters
+ `z`  physical coordinate

## Returns
+ `values`  streamfunction values

## Discussion

  The returned array has one row per `z` value and one column
  per retained wavenumber in `k`.
