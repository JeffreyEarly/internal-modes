---
layout: default
title: modesAtFrequency
parent: InternalModes
grand_parent: Classes
nav_order: 20
mathjax: true
---

#  modesAtFrequency

Return vertical modes for a fixed frequency.


---

## Declaration
```matlab
 [F,G,h,k,varargout] = modesAtFrequency(self,omega,varargin)
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


