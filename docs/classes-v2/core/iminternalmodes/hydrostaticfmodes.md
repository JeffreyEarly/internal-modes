---
layout: default
title: hydrostaticFModes
parent: IMInternalModes
grand_parent: Core
nav_order: 10
mathjax: true
---

#  hydrostaticFModes

Create the hydrostatic `F` internal-mode EVP.


---

## Declaration
```matlab
 evp = IMInternalModes.hydrostaticFModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface boundary condition
+ `options.bottomBoundary`  bottom boundary condition

## Returns
+ `evp`  hydrostatic `F` EVP

## Discussion

  This factory creates the hydrostatic `F`-form problem
  $$-\frac{\partial}{\partial z}
  \left(
  \frac{1}{N^2(z)}
  \frac{\partial F_j}{\partial z}(z)
  \right)
  =
  \lambda_j\frac{F_j(z)}{g},
  \qquad \lambda_j=\frac{1}{h_j}.$$
  At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
  `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
  $$-\left[
  a_\ell F_j(z_\ell)
  -
  b_\ell\frac{1}{N^2(z_\ell)}
  \frac{\partial F_j}{\partial z}(z_\ell)
  \right]
  =
  \lambda_j\left[
  c_\ell F_j(z_\ell)
  -
  d_\ell\frac{1}{N^2(z_\ell)}
  \frac{\partial F_j}{\partial z}(z_\ell)
  \right].$$
  The default surface and bottom boundary conditions are
  `IMBoundaryCondition.neumann()`, giving
  $$\frac{1}{N^2(z_s)}
  \frac{\partial F_j}{\partial z}(z_s)=0,\qquad
  \frac{1}{N^2(z_b)}
  \frac{\partial F_j}{\partial z}(z_b)=0.$$
  Through the hydrostatic relation
  $$G_j(z)=-\frac{g}{N^2(z)}
  \frac{\partial F_j}{\partial z}(z),$$
  these are the same rigid-lid and rigid-bottom conditions
  $$G_j(z_s)=0,\qquad G_j(z_b)=0.$$
  The barotropic zero mode is inferred from the canonical left
  problem during mode selection.
  This factory sets `parameters.formulation` and `parameters.g`;
  `parameters.f0` is supplied by the internal-mode constructor
  default.
