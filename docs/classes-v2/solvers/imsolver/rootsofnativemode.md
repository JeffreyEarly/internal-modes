---
layout: default
title: rootsOfNativeMode
parent: IMSolver
grand_parent: Solvers
nav_order: 14
mathjax: true
---

#  rootsOfNativeMode

Return physical roots of one native mode.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 zRoots = rootsOfNativeMode(solver,nativeMode)
```
## Parameters
+ `nativeMode`  one native mode column

## Returns
+ `zRoots`  physical roots in the solver domain

## Discussion

  Concrete solvers may implement this developer hook when their
  native representation supports accurate root finding. The base
  implementation reports that mode-root grids are unavailable.
