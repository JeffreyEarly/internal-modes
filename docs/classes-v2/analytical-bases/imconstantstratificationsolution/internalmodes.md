---
layout: default
title: internalModes
parent: IMConstantStratificationSolution
grand_parent: Analytical bases
nav_order: 5
mathjax: true
---

#  internalModes

Create an exact internal-mode basis.


---

## Declaration
```matlab
 basisSet = internalModes(solution,evp,options)
```
## Parameters
+ `evp`  internal-mode EVP
+ `options.nModes`  number of retained modes
+ `options.normalization`  active normalization
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  exact analytical internal-mode basis

## Discussion

  Generalized-energy APV modes are returned in ascending
  eigenvalue order $$1/h$$. Hyperbolic negative branches come
  first, followed by an optional affine zero branch with
  `h=Inf`, then trigonometric positive branches. The same
  positive depth-normalization factor scales exact `F` and `G`.
