---
layout: default
title: IMQuadratureWeightFit
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMQuadratureWeightFit

Create diagnostics for a quadrature-weight fit.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 weightFit = IMQuadratureWeightFit(options)
```
## Parameters
+ `options.transform`  transform using fitted weights
+ `options.geometricTransform`  reference transform using geometric weights
+ `options.objectiveName`  least-squares objective name
+ `options.objectiveMatrix`  least-squares matrix
+ `options.objectiveTarget`  least-squares target
+ `options.nonnegativeConstraint`  whether nonnegative weights were imposed
+ `options.depthConstraint`  whether exact depth was imposed
+ `options.depthTarget`  full physical depth
+ `options.exitFlag`  optimizer exit flag
+ `options.solverOutput`  optimizer diagnostics

## Returns
+ `weightFit`  initialized quadrature-weight fit diagnostics

## Discussion

  This developer constructor combines transforms built on the
  same fixed points and retained modes with the least-squares
  system and optimizer result that produced the fitted weights.
  Ordinary users obtain this object from
  `IMBasisSet.quadratureWeightsForPoints`.
