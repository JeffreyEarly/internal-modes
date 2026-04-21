---
layout: default
title: ModesFromGEP
parent: InternalModesFiniteDifference
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  ModesFromGEP

Solve a generalized EVP and map its modes onto the public grid.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [F,G,h,varargout] = ModesFromGEP(self,A,B,h_func,varargin)
```
## Parameters
+ `self`  InternalModesFiniteDifference instance
+ `A`  left generalized-eigenproblem matrix
+ `B`  right generalized-eigenproblem matrix
+ `h_func`  map from eigenvalue to equivalent depth
+ `varargin`  optional requests among `F2`, `G2`, `N2G2`, `uMax`, `wMax`, `kConstant`, and `omegaConstant`

## Returns
+ `F`  horizontal-velocity mode matrix on `zOut`
+ `G`  vertical-velocity mode matrix on `zOut`
+ `h`  equivalent-depth row vector
+ `varargout`  requested diagnostics from `NormalizeModes`

## Discussion


