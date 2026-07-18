---
layout: default
title: fromSolverAtFrequency
parent: InternalModesBasis
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  fromSolverAtFrequency

Solve modes at fixed frequency and return an annotated basis.


---

## Declaration
```matlab
 basis = InternalModesBasis.fromSolverAtFrequency(solver,omega,options)
```
## Parameters
+ `solver`  InternalModesBase solver instance
+ `omega`  fixed frequency in radians per second
+ `options.nModes`  number of mode columns retained
+ `options.useModeAdaptedGrid`  true to request mode-adapted quadrature points when available
+ `options.nQuadraturePoints`  number of quadrature points for mode-adapted solves
+ `options.g`  gravitational acceleration used when solver state does not expose it publicly

## Returns
+ `basis`  InternalModesBasis containing solved modes and component role

## Discussion

For `omega=0`, the returned basis is marked geostrophic and
both F and G components are canonical vertical projection
bases.
