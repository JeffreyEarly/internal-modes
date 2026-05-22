---
layout: default
title: waveModesAtFrequency
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 23
mathjax: true
---

#  waveModesAtFrequency

Create the wave-mode `G` EVP at fixed frequency.


---

## Declaration
```matlab
 evp = IMEigenvalueProblem.waveModesAtFrequency(options)
```
## Parameters
+ `options.omega`  fixed frequency in radians per second
+ `options.f0`  Coriolis parameter in radians per second
+ `options.g`  gravitational acceleration in meters per second squared
+ `options.surfaceBoundary`  location-free surface boundary law
+ `options.bottomBoundary`  location-free bottom boundary law

## Returns
+ `evp`  fixed-frequency wave-mode `G` EVP

## Discussion

  This factory fixes the frequency `omega=options.omega` and
  solves the physical-coordinate strong form
  $$G_{zz}=\lambda(\omega^2-N^2)G/g,\qquad h=1/\lambda.$$
  The solved variable is `G`; the linked diagnostic variable is
  $$F=hG_z.$$
  The factory stores `parameters.omega`, uses the default
  `Normalization.omegaConstant` normalization, and places the
  supplied location-free boundary laws on the surface and bottom.
  Omitted boundaries are rigid `G=0` boundaries. The `kConstant`
  normalization remains available for fixed-frequency basis sets
  and includes boundary trace terms when the boundary laws provide
  them.

  ```matlab
  solver = IMSolverSpectral(N2=@(z) 1e-5*ones(size(z)), zDomain=[-500 0], nEVP=48);
  evp = IMEigenvalueProblem.waveModesAtFrequency(omega=1.2e-3, f0=1e-4);
  basisSet = solver.solveEVP(evp, nModes=6);
  F = basisSet.F(linspace(-500, 0, 100).');
  ```
