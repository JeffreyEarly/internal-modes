---
layout: default
title: partialWindowModes
parent: IMBasisSet
grand_parent: Core
nav_order: 14
mathjax: true
---

#  partialWindowModes

Diagonalize a partial scalar Gram matrix.


---

## Declaration
```matlab
 windowModes = partialWindowModes(basisSet,zMin,zMax)
```
## Parameters
+ `zMin`  lower physical bound
+ `zMax`  upper physical bound

## Returns
+ `windowModes`  window-mode decomposition

## Discussion

  This computes the eigendecomposition of the symmetric
  partial-domain Gram matrix and sorts window-mode eigenvalues
  from largest to smallest.
