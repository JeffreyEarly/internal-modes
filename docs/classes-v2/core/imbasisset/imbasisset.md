---
layout: default
title: IMBasisSet
parent: IMBasisSet
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  IMBasisSet

Create a solved basis set.


---

## Declaration
```matlab
 basisSet = IMBasisSet(options)
```
## Parameters
+ `options.solver`  solver reference
+ `options.evp`  EVP descriptor
+ `options.nativeModes`  native mode columns
+ `options.eigenvalues`  retained eigenvalues
+ `options.h`  equivalent depths
+ `options.modeNumber`  physical mode numbers
+ `options.index`  index summary
+ `options.normalization`  active normalization; omitted uses the EVP default
+ `options.metadata`  additional metadata
+ `options.zDomain`  physical vertical domain
+ `options.N2Function`  buoyancy frequency squared function

## Returns
+ `basisSet`  initialized basis set

## Discussion
