---
layout: default
title: IMAnalyticalInternalModesBasis
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 3
mathjax: true
---

#  IMAnalyticalInternalModesBasis

Create an exact internal-mode basis.


---

## Declaration
```matlab
 basisSet = IMAnalyticalInternalModesBasis(options)
```
## Parameters
+ `options.solution`  analytical solution family
+ `options.evp`  internal-mode EVP
+ `options.h`  equivalent depths
+ `options.modeNumber`  retained-mode labels
+ `options.N2`  buoyancy frequency squared function
+ `options.rawVariableFunction`  exact `F`/`G` evaluator
+ `options.rawUzFunction`  exact solved-variable derivative evaluator
+ `options.normalization`  active normalization
+ `options.metadata`  creation metadata

## Returns
+ `basisSet`  exact internal-mode basis

## Discussion
