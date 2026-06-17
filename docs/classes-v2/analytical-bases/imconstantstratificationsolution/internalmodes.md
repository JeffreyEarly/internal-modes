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
