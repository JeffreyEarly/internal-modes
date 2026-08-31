---
layout: default
title: IMMeanDensityAnomalyModesBasis
parent: IMMeanDensityAnomalyModesBasis
grand_parent: Core
nav_order: 2
mathjax: true
---

#  IMMeanDensityAnomalyModesBasis

Create an aligned mean-density-anomaly basis set.


---

## Declaration
```matlab
 basisSet = IMMeanDensityAnomalyModesBasis(options)
```
## Parameters
+ `options.solver`  configured solver
+ `options.evp`  mean-density-anomaly EVP descriptor
+ `options.nativeModes`  native solved `G` columns
+ `options.eigenvalues`  retained eigenvalues
+ `options.modeNumber`  physical mode labels
+ `options.modeSelectionDiagnostics`  selection diagnostics
+ `options.metadata`  captured problem metadata

## Returns
+ `basisSet`  aligned mean-density-anomaly basis set

## Discussion
