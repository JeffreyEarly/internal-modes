---
layout: default
title: partialWindowModes
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 13
mathjax: true
---

#  partialWindowModes

Diagonalize a partial-depth Gram matrix for `F` or `G`.


---

## Declaration
```matlab
 windowModes = partialWindowModes(basisSet,options)
```
## Parameters
+ `options.variable`  `"F"` or `"G"`
+ `options.zBounds`  integration bounds `[zMin zMax]`

## Returns
+ `windowModes`  window-mode decomposition

## Discussion

This computes the eigendecomposition of the symmetric Gram
matrix on `zBounds`. If `variable` is omitted, the solved
formulation is used. The requested variable must have a known
inner product.
