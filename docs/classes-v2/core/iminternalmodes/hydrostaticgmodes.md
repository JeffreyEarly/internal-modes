---
layout: default
title: hydrostaticGModes
parent: IMInternalModes
grand_parent: Core
nav_order: 11
mathjax: true
---

#  hydrostaticGModes

Create the hydrostatic `G` internal-mode EVP.


---

## Declaration
```matlab
 evp = IMInternalModes.hydrostaticGModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface boundary condition
+ `options.bottomBoundary`  bottom boundary condition

## Returns
+ `evp`  hydrostatic `G` EVP

## Discussion

  This factory creates the hydrostatic `G`-form problem
  $$-\frac{\partial^2 G_j}{\partial z^2}(z)
  =
  \lambda_j\frac{N^2(z)}{g}G_j(z),
  \qquad \lambda_j=\frac{1}{h_j}.$$
  At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
  `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
  $$-\left[
  a_\ell G_j(z_\ell)
  -
  b_\ell\frac{\partial G_j}{\partial z}(z_\ell)
  \right]
  =
  \lambda_j\left[
  c_\ell G_j(z_\ell)
  -
  d_\ell\frac{\partial G_j}{\partial z}(z_\ell)
  \right].$$
  The default surface and bottom boundary conditions are
  `IMBoundaryCondition.dirichlet()`, giving rigid-lid and
  rigid-bottom conditions
  $$G_j(z_s)=0,\qquad G_j(z_b)=0.$$
  Solved hydrostatic basis sets install the `geostrophic`
  normalization rule and use it by default because they set
  `modeFamily` to `"geostrophic"`. This factory sets
  `parameters.formulation`, `parameters.f0`, and `parameters.g`.

  ```matlab
  evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
  solver = IMSolverSpectral(nEVP=128);
  basisSet = solver.solveEVP(evp,nModes=4);
  G = basisSet.G(z);
  ```
