---
layout: default
title: psiz
parent: IMSurfaceGeostrophicModesBasis
grand_parent: Core
nav_order: 9
mathjax: true
---

#  psiz

Evaluate vertical derivatives of SQG streamfunction modes.


---

## Declaration
```matlab
 values = psiz(basisSet,z)
```
## Parameters
+ `z`  physical coordinate

## Returns
+ `values`  streamfunction derivative values

## Discussion

  The returned array has one row per `z` value and one column
  per retained wavenumber in `k`.
