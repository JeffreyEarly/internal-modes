---
layout: default
title: IMBasisSet
parent: IMBasisSet
grand_parent: Core
nav_order: 1
mathjax: true
---

#  IMBasisSet

Create a solved scalar basis set.


---

## Declaration
```matlab
 basisSet = IMBasisSet(options)
```
## Parameters
+ `options.solver`  solver reference
+ `options.evp`  canonical EVP descriptor
+ `options.nativeModes`  native mode columns
+ `options.eigenvalues`  retained eigenvalues
+ `options.modeNumber`  retained-mode labels
+ `options.modeSelectionDiagnostics`  mode-selection diagnostics
+ `options.normalization`  active normalization rule name
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  solved scalar basis set

## Discussion
