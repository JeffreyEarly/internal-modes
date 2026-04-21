---
layout: default
title: FiniteDifferenceMatrix
parent: InternalModesFiniteDifference
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  FiniteDifferenceMatrix

Build a finite-difference differentiation matrix with boundary stencils.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 D = FiniteDifferenceMatrix(numDerivs,x,leftBCDerivs,rightBCDerivs,orderOfAccuracy)
```
## Parameters
+ `numDerivs`  derivative order to approximate
+ `x`  grid on which to build the stencil matrix
+ `leftBCDerivs`  derivative order enforced at the left boundary
+ `rightBCDerivs`  derivative order enforced at the right boundary
+ `orderOfAccuracy`  formal order of accuracy

## Returns
+ `D`  differentiation matrix

## Discussion


