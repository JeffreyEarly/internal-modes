---
layout: default
title: coordinateProfile
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 5
mathjax: true
---

#  coordinateProfile

Return fields needed by a solver coordinate map.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 profile = coordinateProfile(evp,coordinateKind)
```
## Parameters
+ `coordinateKind`  solver coordinate kind

## Returns
+ `profile`  struct with coordinate resources

## Discussion

  Generic canonical EVPs are independent of stratification and
  only support the physical `z` coordinate.
