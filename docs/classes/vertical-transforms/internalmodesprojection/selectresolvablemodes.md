---
layout: default
title: selectResolvableModes
parent: InternalModesProjection
grand_parent: Classes
nav_order: 26
mathjax: true
---

#  selectResolvableModes

Select observation-resolvable columns with pivoted QR.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [retainedModes,diagnostics] = InternalModesProjection.selectResolvableModes(B,weights,rankTolerance,maxConditionNumber)
```
## Parameters
+ `B`  sampled observation matrix
+ `weights`  observation weights
+ `rankTolerance`  relative QR pivot tolerance
+ `maxConditionNumber`  maximum retained Gram condition number

## Returns
+ `retainedModes`  selected column numbers
+ `diagnostics`  structure with QR diagnostics

## Discussion
