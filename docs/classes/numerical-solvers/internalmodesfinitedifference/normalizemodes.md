---
layout: default
title: NormalizeModes
parent: InternalModesFiniteDifference
grand_parent: Classes
nav_order: 11
mathjax: true
---

#  NormalizeModes

Normalize finite-difference modes using the active convention.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [F,G,varargout] = NormalizeModes(self,F,G,N2,z,varargin)
```
## Parameters
+ `self`  InternalModesFiniteDifference instance
+ `F`  horizontal-velocity mode matrix
+ `G`  vertical-velocity mode matrix
+ `N2`  buoyancy-frequency profile associated with the modes
+ `z`  depth grid associated with the modes
+ `varargin`  optional requests among `F2`, `G2`, `N2G2`, `uMax`, `wMax`, `kConstant`, and `omegaConstant`

## Returns
+ `F`  normalized horizontal-velocity mode matrix
+ `G`  normalized vertical-velocity mode matrix
+ `varargout`  requested quadratic-integral and normalization diagnostics

## Discussion


