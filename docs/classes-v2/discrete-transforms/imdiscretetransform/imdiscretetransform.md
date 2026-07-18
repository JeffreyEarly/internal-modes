---
layout: default
title: IMDiscreteTransform
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMDiscreteTransform

Create a scalar discrete transform from canonical matrices.


---

## Declaration
```matlab
 transform = IMDiscreteTransform(options)
```
## Parameters
+ `options.z`  physical sample points
+ `options.increments`  quadrature increments
+ `options.modeNumber`  retained mode labels
+ `options.normalization`  basis normalization name
+ `options.inverseMatrix`  inverse transform matrix containing the sampled modes
+ `options.metricMatrix`  sampled metric matrix
+ `options.targetGramMatrix`  continuous diagonal Gram target

## Returns
+ `transform`  initialized scalar discrete transform

## Discussion

Ordinary users construct this object with
`IMBasisSet.discreteTransform`. This constructor is the
canonical matrix-level initialization path for later discrete
transform builders.
