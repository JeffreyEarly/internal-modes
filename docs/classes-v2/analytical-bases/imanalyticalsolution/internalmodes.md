---
layout: default
title: internalModes
parent: IMAnalyticalSolution
grand_parent: Analytical bases
nav_order: 6
mathjax: true
---

#  internalModes

Create an exact internal-mode basis when available.


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

  Concrete solution families return an
  `IMAnalyticalInternalModesBasis` for supported EVPs and throw a
  class-specific unavailable or unsupported error otherwise.
