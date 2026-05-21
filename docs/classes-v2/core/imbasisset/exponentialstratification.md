---
layout: default
title: exponentialStratification
parent: IMBasisSet
grand_parent: Classes
nav_order: 13
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
+ `options.evp`  supported rigid-endpoint `G` eigenvalue-problem descriptor
+ `options.N0`  surface buoyancy frequency
+ `options.b`  exponential e-folding depth
+ `options.zDomain`  physical vertical domain with surface at zero
+ `options.nModes`  number of modes
+ `options.normalization`  active normalization; omitted uses the EVP default
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  analytical exponential-stratification basis set

## Discussion

  The returned basis set evaluates exact rigid-endpoint
  `G`-formulation modes for
  $$N^2(z)=N_0^2e^{2z/b}$$ without solving a numerical EVP.
