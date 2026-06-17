---
layout: default
title: IMInternalModesBasis
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 3
mathjax: true
---

#  IMInternalModesBasis

Create an internal-mode basis set.


---

## Declaration
```matlab
 basisSet = IMInternalModesBasis(options)
```
## Parameters
+ `options.solver`  solver reference
+ `options.evp`  internal-mode EVP descriptor
+ `options.nativeModes`  native mode columns
+ `options.eigenvalues`  retained eigenvalues
+ `options.h`  equivalent depths
+ `options.modeNumber`  physical mode numbers
+ `options.modeSelectionDiagnostics`  mode-selection diagnostics
+ `options.normalization`  active normalization rule name or enum value
+ `options.metadata`  additional metadata
+ `options.zDomain`  physical vertical domain

## Returns
+ `basisSet`  internal-mode basis set

## Discussion
