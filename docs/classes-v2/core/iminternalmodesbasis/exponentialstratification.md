---
layout: default
title: exponentialStratification
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 6
mathjax: true
---

#  exponentialStratification

Create an analytical exponential-stratification basis set.


---

## Declaration
```matlab
 basisSet = IMInternalModesBasis.exponentialStratification(options)
```
## Parameters
+ `options.evp`  internal-mode EVP descriptor
+ `options.N0`  surface buoyancy frequency
+ `options.b`  exponential e-folding depth
+ `options.zDomain`  physical vertical domain
+ `options.nModes`  number of modes
+ `options.normalization`  active normalization
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  analytical exponential-stratification basis set

## Discussion
