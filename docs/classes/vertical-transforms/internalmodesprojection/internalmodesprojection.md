---
layout: default
title: InternalModesProjection
parent: InternalModesProjection
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  InternalModesProjection

Create an observation projection from canonical persisted state.


---

## Declaration
```matlab
 projection = InternalModesProjection(options)
```
## Parameters
+ `options.component`  `"F"` or `"G"`
+ `options.projectionMatrix`  matrix from observations to retained coefficients
+ `options.resolutionMatrix`  matrix from candidate coefficients to recovered coefficients
+ `options.spectralWindow`  spectrum transfer matrix

## Returns
+ `projection`  initialized InternalModesProjection instance

## Discussion
