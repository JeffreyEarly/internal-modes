---
layout: default
title: innerProduct
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 16
mathjax: true
---

#  innerProduct

Return the declared inner-product recipe for a variable.


---

## Declaration
```matlab
 spec = innerProduct(evp,variable)
```
## Parameters
+ `variable`  `"F"` or `"G"`

## Returns
+ `spec`  struct with `variable`, `interiorWeight`, `surfaceWeights`, `bottomWeights`, and `hasKnownBoundaryWeights`

## Discussion

  The returned struct contains the interior weight and the
  endpoint weights whose `innerProduct` matches `variable`. The
  endpoint factors inside those weights may still evaluate either
  `F` or `G`; the selector is the inner product being constructed.
