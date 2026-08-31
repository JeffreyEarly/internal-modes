---
layout: default
title: fixedGridTransform
parent: InternalModesBasis
grand_parent: Classes
nav_order: 10
mathjax: true
---

#  fixedGridTransform

Build a transform on another transform's vertical grid.


---

## Declaration
```matlab
 transform = fixedGridTransform(self,referenceTransform,options)
```
## Parameters
+ `self`  InternalModesBasis instance
+ `referenceTransform`  InternalModesTransform whose grid and weights define the fixed grid
+ `options.component`  `"F"`, `"G"`, or `"both"`
+ `options.projectionTolerance`  leakage tolerance for retained prefix modes
+ `options.maxConditionNumber`  condition-number limit
+ `options.nTailCheck`  number of rejected tail modes checked for leakage
+ `options.preserveSize`  true to keep rejected rows as zeros

## Returns
+ `transform`  InternalModesTransform evaluated on the fixed grid

## Discussion

  This method evaluates the current basis on the reference grid
  and builds a weighted pseudoinverse there. It is intended for
  cases where a model fixes the hydrostatic quadrature grid but
  needs nonzero-$$\kappa$$ IGW modes evaluated on that grid.
