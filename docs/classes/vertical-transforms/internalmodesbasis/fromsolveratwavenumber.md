---
layout: default
title: fromSolverAtWavenumber
parent: InternalModesBasis
grand_parent: Classes
nav_order: 15
mathjax: true
---

#  fromSolverAtWavenumber

Solve modes at fixed horizontal wavenumber.


---

## Declaration
```matlab
 basis = InternalModesBasis.fromSolverAtWavenumber(solver,kappa,options)
```
## Parameters
+ `solver`  InternalModesBase solver instance
+ `kappa`  horizontal wavenumber in radians per meter
+ `options.nModes`  number of mode columns retained
+ `options.useModeAdaptedGrid`  true to request mode-adapted quadrature points when available
+ `options.nQuadraturePoints`  number of quadrature points for mode-adapted solves
+ `options.g`  gravitational acceleration used when solver state does not expose it publicly

## Returns
+ `basis`  InternalModesBasis containing solved modes and component role

## Discussion

For nonzero $$\kappa$$, the returned basis marks G as the
canonical Sturm-Liouville projection component and F as
diagnostic/evaluation-only.
