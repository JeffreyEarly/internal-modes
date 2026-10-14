---
layout: default
title: IMGeostrophicModesBasis
parent: IMGeostrophicModesBasis
grand_parent: Core
nav_order: 2
mathjax: true
---

#  IMGeostrophicModesBasis

Create a basis from aligned native `F` and `G` columns.


---

## Declaration
```matlab
 basisSet = IMGeostrophicModesBasis(options)
```
## Parameters
+ `options.solver`  configured solver
+ `options.evp`  geostrophic EVP descriptor
+ `options.nativeModes`  native solved `F` columns
+ `options.nativeGModes`  aligned native solved `G` columns
+ `options.companionSolver`  optional G discretization for staggered finite differences
+ `options.eigenvalues`  retained public eigenvalues
+ `options.inverseDepth`  retained solved inverse depths
+ `options.modeNumber`  physical mode labels
+ `options.modeSelectionDiagnostics`  selection diagnostics
+ `options.normalization`  normalization rule name or enum value
+ `options.metadata`  captured problem metadata

## Returns
+ `basisSet`  aligned geostrophic basis set

## Discussion
