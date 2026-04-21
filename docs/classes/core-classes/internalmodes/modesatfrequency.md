---
layout: default
title: ModesAtFrequency
parent: InternalModes
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  ModesAtFrequency

Return vertical modes for a fixed frequency.


---

## Declaration
```matlab
 [F,G,h,k,varargout] = ModesAtFrequency(self,omega,varargin)
```
## Parameters
+ `self`  InternalModes instance
+ `omega`  frequency in radians per second
+ `varargin`  additional requests forwarded to the concrete solver

## Returns
+ `F`  horizontal-velocity mode matrix on `zOut`
+ `G`  vertical-velocity mode matrix on `zOut`
+ `h`  equivalent-depth row vector
+ `k`  horizontal wavenumber row vector implied by `h` and `omega`
+ `varargout`  additional outputs forwarded from the concrete solver

## Discussion


