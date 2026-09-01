---
layout: default
title: targetMajorantGramMatrix
parent: IMInternalModesDiscreteTransform
grand_parent: Discrete transforms
nav_order: 30
mathjax: true
---

#  targetMajorantGramMatrix

Return the continuous positive Hilbert-majorant Gram matrix.


---

## Declaration
```matlab
 matrix = targetMajorantGramMatrix(transform,options)
```
## Parameters
+ `options.variable`  `"F"` or `"G"`

## Returns
+ `matrix`  continuous positive majorant Gram matrix

## Discussion

  `targetGramMatrix` remains the signed Pontryagin pairing used
  for projection. This accessor returns the positive matrix used
  for error magnitudes and coupled quadratic certification.
