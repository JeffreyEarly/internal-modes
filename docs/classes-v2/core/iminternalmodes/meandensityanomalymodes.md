---
layout: default
title: meanDensityAnomalyModes
parent: IMInternalModes
grand_parent: Core
nav_order: 16
mathjax: true
---

#  meanDensityAnomalyModes

Create generalized-energy mean-density-anomaly modes.


---

## Declaration
```matlab
 evp = IMInternalModes.meanDensityAnomalyModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.g`  gravitational acceleration
+ `options.g0`  signed finite surface acceleration, zero, or positive infinity
+ `options.gd`  signed finite bottom acceleration, zero, or positive infinity

## Returns
+ `evp`  generalized-energy mean-density-anomaly EVP

## Discussion

  This factory creates the `G`-form problem

  $$
  -G_j''(z)=\frac{N^2(z)}{g h_j}G_j(z)
  $$

  with endpoint conditions

  $$
  g h_jG_j'(z_s)=g_0G_j(z_s),\qquad
  g h_jG_j'(z_b)=-g_dG_j(z_b).
  $$

  `g0` and `gd` are required. A finite value, including zero,
  keeps that endpoint active; zero is the Neumann limit.
  Positive infinity imposes Dirichlet data and omits the
  corresponding generalized-energy endpoint term. `NaN` and
  negative infinity are rejected.

  Solved modes use the signed generalized-energy normalization

  $$
  \frac{1}{g}\int_{z_b}^{z_s}N^2G_iG_j\,dz
  +\frac{g_0}{g}G_i(z_s)G_j(z_s)
  +\frac{g_d}{g}G_i(z_b)G_j(z_b)
  =\epsilon_j\delta_{ij},
  $$

  with inactive terms omitted. The basis exposes
  `basisSet.signatures` as $$\epsilon_j\in\{-1,+1\}$$.
  For the continuous projection functional
  $$\mathcal G_j[X]=\langle G_j,X\rangle_G$$, normalized
  coefficients obey $$A_j=\epsilon_j\mathcal G_j[X]$$.
  Its aligned diagnostic pressure modes are computed by
  surface-referenced integration,

  $$
  F_j(z)=\frac{1}{g}\int_z^{z_s}N^2(z')G_j(z')\,dz',
  \qquad F_j(z_s)=0.
  $$

  Discrete transforms directly project the `G` channel and
  synthesize both `G` and `F`. The diagnostic `F` channel does
  not define an independent coefficient projection metric.

  ```matlab
  evp = IMInternalModes.meanDensityAnomalyModes( ...
      N2=N2,zDomain=[-4000 0],g0=0.02,gd=Inf);
  basisSet = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=8);
  F = basisSet.F(z);
  G = basisSet.G(z);
  ```
