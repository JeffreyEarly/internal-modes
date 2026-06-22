---
layout: default
title: F
parent: IMSurfaceGeostrophicModesBasis
grand_parent: Core
nav_order: 1
mathjax: true
---

#  F

Evaluate projected SQG streamfunction modes.


---

## Declaration
```matlab
 values = F(basisSet,z)
```
## Parameters
+ `z`  physical coordinate

## Returns
+ `values`  projected `F` values

## Discussion

  The returned array has one row per `z` value and one column
  per projected boundary mode.
