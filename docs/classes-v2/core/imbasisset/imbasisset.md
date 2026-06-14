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
+ `options.h`  equivalent depths
+ `options.modeNumber`  physical mode numbers
+ `options.index`  selection diagnostics
+ `options.normalization`  active normalization
+ `options.metadata`  additional metadata
+ `options.zDomain`  physical vertical domain
+ `options.N2Function`  buoyancy frequency squared function

## Returns
+ `basisSet`  solved scalar basis set

## Discussion
