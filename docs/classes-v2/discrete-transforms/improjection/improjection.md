---
layout: default
title: IMProjection
parent: IMProjection
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMProjection

Prepare a projection without making acceptance decisions.


---

## Declaration
```matlab
 projection = IMProjection(sampledBasis,metricMatrix,targetGramMatrix,options)
```
## Parameters
+ `sampledBasis`  real finite sample-by-column matrix
+ `metricMatrix`  real symmetric sample pairing
+ `targetGramMatrix`  real symmetric continuous pairing
+ `options.majorantGramMatrix`  positive active coefficient metric
+ `options.activeColumnMask`  logical row, default all columns
+ `options.columnLabels`  string row of scientific labels; default positional labels
+ `options.provenance`  construction identity and numerical evidence

## Returns
+ `self`  fixed projection and numerical diagnostics

## Discussion

  A diagonal target defaults to its absolute diagonal majorant.
  A general target requires an explicit positive majorant. Active
  target diagonal entries must be nonzero for Gram normalization.
  Inactive columns must be numerically zero in sampled and target
  data; their forward rows remain zero. Empty active sets are valid.
  With no arguments, construct the canonical empty Galerkin projection.
