---
layout: default
title: constantStratification
parent: IMBasisSet
grand_parent: Classes
nav_order: 5
mathjax: true
---

#  constantStratification

Create an analytical constant-stratification basis set.


---

## Declaration
```matlab
 basisSet = IMBasisSet.constantStratification(options)
```
## Parameters
+ `options.evp`  supported wave-mode or hydrostatic eigenvalue-problem descriptor
+ `options.N0`  constant buoyancy frequency
+ `options.zDomain`  physical vertical domain
+ `options.nModes`  number of modes
+ `options.normalization`  active normalization; omitted uses the EVP default
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  analytical constant-stratification basis set

## Discussion

  The returned basis set evaluates exact constant-stratification
  `F` and `G` modes without solving a numerical EVP.
