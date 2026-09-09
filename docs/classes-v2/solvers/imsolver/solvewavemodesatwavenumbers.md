---
layout: default
title: solveWaveModesAtWavenumbers
parent: IMSolver
grand_parent: Solvers
nav_order: 19
mathjax: true
---

#  solveWaveModesAtWavenumbers

Solve continuous wave bases for an ordered collection of wavenumbers.


---

## Declaration
```matlab
 [collection,diagnostics] = solveWaveModesAtWavenumbers(solver,kappa,options)
```
## Parameters
+ `kappa`  nonempty row of requested nonnegative wavenumbers, in radians per meter
+ `options.N2`  buoyancy frequency squared evaluator
+ `options.zDomain`  increasing physical vertical bounds
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface condition of the G EVP
+ `options.bottomBoundary`  bottom condition of the G EVP
+ `options.nModes`  scalar or row of positive counts aligned with kappa(kappa>0); empty only for all-zero kappa
+ `options.nInertialModes`  independent count, required when zero is requested

## Returns
+ `collection`  continuous bases and exact requested-page mapping
+ `diagnostics`  construction costs and solve provenance, not accuracy qualification

## Discussion

  Each distinct requested wavenumber is solved exactly once. Coordinates,
  stratification samples, derivative matrices and the shared pencil terms
  are prepared once for this call. The wave pencil is
  $$A(\kappa)=A_0+\kappa^2 A_2,\qquad B(\kappa)=B_0,$$
  with zero endpoint rows in $$A_2$$ because boundary equations replace the
  interior operator there. No approximate sharing or global cache is used.
  Shared preparation supports the exact built-in IMSolverSpectral and
  IMSolverFiniteDifference classes. Custom solvers must use independent
  solveEVP calls until they expose a validated preparation-reuse capability.

  Positive wavenumbers request `nModes` columns, either a uniform scalar or
  one count per positive entry of `kappa`, preserving its order and repeats.
  Repeated wavenumbers must request identical counts. A zero-wavenumber request
  requires an explicit independent `nInertialModes` count. Both counts are
  exact requests: an insufficient returned family raises an error. Stored
  bases retain the same scientific labels, normalization and boundary
  conditions as independent `solveEVP` calls. Frequency signs and physical
  wave polarizations belong to the caller, not this vertical-basis solve.

  ```matlab
  [bases,costs] = solver.solveWaveModesAtWavenumbers([2e-4 0 1e-4 2e-4],N2=N2,zDomain=[-1000 0],nModes=8,nInertialModes=5);
  G = bases.evaluate(z,variable="G",pages=[1 3]);
  ```
