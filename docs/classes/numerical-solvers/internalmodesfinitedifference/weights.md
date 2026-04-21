---
layout: default
title: weights
parent: InternalModesFiniteDifference
grand_parent: Classes
nav_order: 17
mathjax: true
---

#  weights

Return Fornberg finite-difference weights for one stencil location.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 c = weights(z,x,m)
```
## Parameters
+ `z`  evaluation location
+ `x`  stencil grid points
+ `m`  highest derivative order requested

## Returns
+ `c`  weight matrix whose rows correspond to derivative order

## Discussion


