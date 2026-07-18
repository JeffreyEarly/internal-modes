---
layout: default
title: nativeTransform
parent: InternalModesBasis
grand_parent: Classes
nav_order: 22
mathjax: true
---

#  nativeTransform

Build a transform on the basis' native vertical grid.


---

## Declaration
```matlab
 transform = nativeTransform(self,options)
```
## Parameters
+ `self`  InternalModesBasis instance
+ `options.component`  `"F"`, `"G"`, or `"both"`
+ `options.nModes`  maximum number of modes retained
+ `options.projectionMethod`  `"auto"`, `"directInverse"`, `"weightedPseudoinverse"`, or `"canonical"`
+ `options.allowNoncanonical`  true to allow numerical wave-F projections
+ `options.maxConditionNumber`  condition-number limit for direct inverses

## Returns
+ `transform`  InternalModesTransform on the native grid

## Discussion

The native transform maps between modal coefficients and
fields sampled at `basis.z`. Components with Dirichlet rows,
such as rigid-lid G modes, are solved on their active interior
rows and expanded back to full-grid matrices. With
`projectionMethod="auto"`, a square, well-conditioned active
inverse matrix uses a direct inverse; otherwise the transform
uses a weighted pseudoinverse.
