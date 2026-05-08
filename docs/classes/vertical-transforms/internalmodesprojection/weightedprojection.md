---
layout: default
title: weightedProjection
parent: InternalModesProjection
grand_parent: Classes
nav_order: 31
mathjax: true
---

#  weightedProjection

Build a weighted projection for selected observation columns.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [A,resolutionMatrix,gramMatrix] = InternalModesProjection.weightedProjection(B,weights,retainedModes)
```
## Parameters
+ `B`  sampled observation matrix
+ `weights`  observation weights
+ `retainedModes`  retained column numbers

## Returns
+ `A`  weighted projection matrix
+ `resolutionMatrix`  coefficient resolution matrix
+ `gramMatrix`  retained weighted Gram matrix

## Discussion
