---
layout: default
title: modesAtWavenumber
parent: InternalModes
grand_parent: Classes
nav_order: 21
mathjax: true
---

#  modesAtWavenumber

Return vertical modes for a fixed horizontal wavenumber.


---

## Declaration
```matlab
 [F,G,h,omega,varargout] = modesAtWavenumber(self,k,varargin)
```
## Parameters
+ `self`  InternalModes instance
+ `k`  horizontal wavenumber
+ `varargin`  additional requests forwarded to the concrete solver

## Returns
+ `F`  horizontal-velocity mode matrix on `zOut`
+ `G`  vertical-velocity mode matrix on `zOut`
+ `h`  equivalent-depth row vector
+ `omega`  frequency row vector implied by `h` and `k`
+ `varargout`  additional outputs forwarded from the concrete solver

## Discussion


