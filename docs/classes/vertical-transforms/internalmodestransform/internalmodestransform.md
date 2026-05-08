---
layout: default
title: InternalModesTransform
parent: InternalModesTransform
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  InternalModesTransform

Create a vertical transform from canonical persisted state.


---

## Declaration
```matlab
 transform = InternalModesTransform(options)
```
## Parameters
+ `options.z`  depth grid
+ `options.forwardF`  F forward matrix
+ `options.inverseF`  F inverse matrix
+ `options.forwardG`  G forward matrix
+ `options.inverseG`  G inverse matrix
+ `options.spectralWeightsF`  F spectrum weights
+ `options.spectralWeightsG`  G spectrum weights

## Returns
+ `transform`  initialized InternalModesTransform instance

## Discussion

  This constructor stores already-built operators. Use
  `InternalModesBasis` factories to solve modes and build new
  transforms from stratification.
