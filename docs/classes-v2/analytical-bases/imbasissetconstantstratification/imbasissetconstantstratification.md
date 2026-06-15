---
layout: default
title: IMBasisSetConstantStratification
parent: IMBasisSetConstantStratification
grand_parent: Analytical bases
nav_order: 1
mathjax: true
---

#  IMBasisSetConstantStratification

Create an exact constant-stratification basis set.


---

## Declaration
```matlab
 basisSet = IMBasisSetConstantStratification(options)
```
## Parameters
+ `options.evp`  supported wave-mode or hydrostatic eigenvalue-problem descriptor
+ `options.N0`  constant buoyancy frequency
+ `options.zDomain`  physical vertical domain
+ `options.nModes`  number of retained modes
+ `options.normalization`  active normalization rule; omitted uses the EVP default
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  exact constant-stratification basis set

## Discussion
