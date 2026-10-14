---
layout: default
title: IMGeostrophicZeroAPVModesBasis
parent: IMGeostrophicZeroAPVModesBasis
grand_parent: Core
nav_order: 3
mathjax: true
---

#  IMGeostrophicZeroAPVModesBasis

Create a canonical boundary-normalized zero-APV basis.


---

## Declaration
```matlab
 basisSet = IMGeostrophicZeroAPVModesBasis(options)
```
## Parameters
+ `options.problem`  geostrophic zero-APV problem
+ `options.solver`  configured numerical solver
+ `options.nativeModes`  canonical native F structures
+ `options.nativeGModes`  optional aligned native G structures
+ `options.companionSolver`  optional G discretization; requires nativeGModes
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  canonical boundary-normalized basis

## Discussion
