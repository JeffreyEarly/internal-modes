---
layout: default
title: partialWindowModes
parent: IMBasisSet
grand_parent: Core
nav_order: 20
mathjax: true
---

#  partialWindowModes

Diagonalize a partial scalar Gram matrix.


---

## Declaration
```matlab
 windowModes = partialWindowModes(basisSet,options)
```
## Parameters
+ `options.zBounds`  integration bounds `[zMin zMax]`

## Returns
+ `windowModes`  window-mode decomposition

## Discussion

This computes the eigendecomposition of the symmetric Gram
matrix on `zBounds` and sorts window-mode eigenvalues from
largest to smallest.
