---
layout: default
title: waveModesAtWavenumber
parent: IMInternalModes
grand_parent: Core
nav_order: 17
mathjax: true
---

#  waveModesAtWavenumber

Create the fixed-wavenumber wave-mode EVP.


---

## Declaration
```matlab
 evp = IMInternalModes.waveModesAtWavenumber(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.k`  horizontal wavenumber
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface boundary condition
+ `options.bottomBoundary`  bottom boundary condition

## Returns
+ `evp`  fixed-wavenumber `G` EVP

## Discussion

This factory creates the fixed-wavenumber `G`-form problem

$$
-\frac{\partial^2 G_j}{\partial z^2}(z)
+k^2G_j(z)
=
\lambda_j\frac{N^2(z)-f_0^2}{g}G_j(z),
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
Solved fixed-wavenumber basis sets install the `kConstant`
normalization rule and use it by default.
This factory adds `parameters.k` and sets
`parameters.formulation`, `parameters.f0`, and `parameters.g`.
