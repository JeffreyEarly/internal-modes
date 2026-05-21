---
layout: default
title: partialDepthPE
parent: IMBoundary
grand_parent: Classes
nav_order: 25
mathjax: true
---

#  partialDepthPE

Create partial-depth potential-energy active boundary conditions.


---

## Declaration
```matlab
 boundaryConditions = IMBoundary.partialDepthPE(options)
```
## Parameters
+ `options.boundarySign`  `"positive"` or `"negative"`

## Returns
+ `boundaryConditions`  bottom and surface active conditions

## Discussion

  Positive boundary signs add no negative index directions.
  Negative boundary signs add one negative direction at each
  window endpoint.
