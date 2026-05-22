---
layout: default
title: IMBasisSetExponentialStratification
parent: IMBasisSetExponentialStratification
grand_parent: Analytical bases
nav_order: 1
mathjax: true
---

#  IMBasisSetExponentialStratification

Create an exact exponential-stratification basis set.


---

## Declaration
```matlab
 basisSet = IMBasisSetExponentialStratification(options)
```
## Parameters
+ `options.evp`  supported rigid-endpoint `G` eigenvalue-problem descriptor
+ `options.N0`  surface buoyancy frequency
+ `options.b`  exponential e-folding depth
+ `options.zDomain`  physical vertical domain with surface at zero
+ `options.nModes`  number of retained modes
+ `options.normalization`  active normalization; omitted uses the EVP default
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  exact exponential-stratification basis set

## Discussion
