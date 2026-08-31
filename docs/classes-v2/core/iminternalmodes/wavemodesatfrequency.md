---
layout: default
title: waveModesAtFrequency
parent: IMInternalModes
grand_parent: Core
nav_order: 18
mathjax: true
---

#  waveModesAtFrequency

Create the fixed-frequency wave-mode EVP.


---

## Declaration
```matlab
 evp = IMInternalModes.waveModesAtFrequency(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.omega`  wave frequency
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface boundary condition
+ `options.bottomBoundary`  bottom boundary condition

## Returns
+ `evp`  fixed-frequency `G` EVP

## Discussion

  This factory creates the fixed-frequency `G`-form problem

  $$
  -\frac{\partial^2 G_j}{\partial z^2}(z)
  =
  \lambda_j\frac{N^2(z)-\omega^2}{g}G_j(z),
  \qquad \lambda_j=\frac{1}{h_j}.
  $$

  At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
  `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as

  $$
  -\left[
  a_\ell G_j(z_\ell)
  -b_\ell\frac{\partial G_j}{\partial z}(z_\ell)
  \right]
  =
  \lambda_j\left[
  c_\ell G_j(z_\ell)
  -d_\ell\frac{\partial G_j}{\partial z}(z_\ell)
  \right].
  $$

  The default surface and bottom boundary conditions are
  `IMBoundaryCondition.dirichlet()`, giving rigid endpoint
  conditions

  $$
  G_j(z_s)=0,\qquad G_j(z_b)=0.
  $$

  A linear free-surface condition at the surface can be written
  as

  $$
  G_j(z_s)=h_j\frac{\partial G_j}{\partial z}(z_s),
  \qquad \lambda_j=\frac{1}{h_j},
  $$

  equivalently

  $$
  \frac{\partial G_j}{\partial z}(z_s)
  =
  \lambda_j G_j(z_s).
  $$

  In canonical boundary-condition coefficients this is
  `IMBoundaryCondition(a=0,b=1,c=1,d=0)` at the surface.
  Solved fixed-frequency basis sets use the generic `unity`
  normalization by default. A fixed-frequency diagnostic `F`
  inner-product normalization is deferred until the wave
  diagnostic inner-product catalog is derived. This factory
  adds `parameters.omega` and sets `parameters.formulation`,
  `parameters.f0`, and `parameters.g`.
