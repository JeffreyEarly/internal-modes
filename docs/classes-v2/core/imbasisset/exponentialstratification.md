---
layout: default
title: exponentialStratification
parent: IMBasisSet
grand_parent: Core
nav_order: 10
mathjax: true
---

#  exponentialStratification

Create an analytical exponential-stratification basis set.


---

## Declaration
```matlab
 basisSet = IMBasisSet.exponentialStratification(options)
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
