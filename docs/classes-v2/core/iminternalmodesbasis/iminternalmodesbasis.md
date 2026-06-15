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
+ `options.index`  selection diagnostics
+ `options.normalization`  active normalization
+ `options.metadata`  additional metadata
+ `options.zDomain`  physical vertical domain
+ `options.N2`  buoyancy frequency squared function

## Returns
+ `basisSet`  internal-mode basis set

## Discussion
