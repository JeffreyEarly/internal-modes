---
layout: default
title: InternalModes
parent: InternalModes
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  InternalModes

Initialize the wrapper from a profile or a built-in benchmark case.


---

## Declaration
```matlab
 im = InternalModes(varargin)
```
## Parameters
+ `varargin`  initialization arguments for either a user-specified profile or a built-in benchmark case

## Returns
+ `im`  wrapper instance that forwards to a concrete solver

## Discussion

  The main calling patterns are:

  - `InternalModes(rho, zIn, zOut, latitude, ...)`
  - `InternalModes(profileName, methodName, nPoints)`

  where additional name-value pairs are split between wrapper
  options such as `method` and the concrete solver
  constructor.


