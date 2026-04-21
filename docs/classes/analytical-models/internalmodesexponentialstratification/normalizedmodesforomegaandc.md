---
layout: default
title: NormalizedModesForOmegaAndC
parent: InternalModesExponentialStratification
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  NormalizedModesForOmegaAndC

Evaluate and normalize analytical mode functions at `zOut`.


---

## Declaration
```matlab
 [F,G,varargout] = NormalizedModesForOmegaAndC(self,omega,c,varargin)
```
## Parameters
+ `self`  InternalModesExponentialStratification instance
+ `omega`  frequency row vector
+ `c`  phase-speed row vector
+ `varargin`  optional requests among `F2`, `G2`, `N2G2`, `uMax`, `wMax`, `kConstant`, and `omegaConstant`

## Returns
+ `F`  normalized horizontal-velocity mode matrix
+ `G`  normalized vertical-velocity mode matrix
+ `varargout`  requested normalization and quadratic-integral diagnostics

## Discussion


