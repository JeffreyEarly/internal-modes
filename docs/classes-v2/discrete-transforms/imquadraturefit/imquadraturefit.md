---
layout: default
title: IMQuadratureFit
parent: IMQuadratureFit
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMQuadratureFit

Create a fixed-point quadrature fit result.


---

## Declaration
```matlab
 fit = IMQuadratureFit(options)
```
## Parameters
+ `options.fittedTransform`  transform using fitted increments
+ `options.geometricTransform`  transform using geometric increments
+ `options.objectiveName`  least-squares objective name
+ `options.objectiveMatrix`  least-squares matrix
+ `options.objectiveTarget`  least-squares target
+ `options.nonnegativeConstraint`  whether nonnegative increments were imposed
+ `options.depthConstraint`  whether exact depth was imposed
+ `options.depthTarget`  full physical depth
+ `options.exitFlag`  optimizer exit flag
+ `options.solverOutput`  optimizer diagnostics

## Returns
+ `fit`  initialized quadrature fit result

## Discussion

Ordinary users construct this object with
`IMBasisSet.fitQuadrature`.
