---
layout: default
title: dzLogN2
parent: IMInternalModes
grand_parent: Core
nav_order: 3
mathjax: true
---

#  dzLogN2

Evaluate the vertical derivative of `log(N2)`.


---

## Declaration
```matlab
 values = dzLogN2(evp,z)
```
## Parameters
+ `z`  physical coordinate

## Returns
+ `values`  derivative values

## Discussion

  This derivative is used by coordinate mappings that need the
  stratification slope. For more than one point it is computed
  by finite differences on the supplied coordinate vector.
