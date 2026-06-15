---
layout: default
title: constantStratification
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 5
mathjax: true
---

#  constantStratification

Create an analytical constant-stratification basis set.


---

## Declaration
```matlab
 basisSet = IMInternalModesBasis.constantStratification(options)
```
## Parameters
+ `options.evp`  internal-mode EVP descriptor
+ `options.N0`  constant buoyancy frequency
+ `options.zDomain`  physical vertical domain
+ `options.nModes`  number of modes
+ `options.normalization`  active normalization
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  analytical constant-stratification basis set

## Discussion
