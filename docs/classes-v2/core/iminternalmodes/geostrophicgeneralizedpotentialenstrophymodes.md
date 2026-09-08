---
layout: default
title: geostrophicGeneralizedPotentialEnstrophyModes
parent: IMInternalModes
grand_parent: Core
nav_order: 10
mathjax: true
---

#  geostrophicGeneralizedPotentialEnstrophyModes

Create free-surface generalized-potential-enstrophy modes.


---

## Declaration
```matlab
 evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function in radians squared per second squared
+ `options.zDomain`  physical vertical domain in meters
+ `options.k`  positive horizontal wavenumber in radians per meter
+ `options.f0`  nonzero Coriolis parameter in radians per second
+ `options.g`  gravitational acceleration in meters per second squared
+ `options.alpha0`  optional positive surface generalized-potential-enstrophy weight
+ `options.alphaD`  positive bottom generalized-potential-enstrophy weight or `Inf`

## Returns
+ `evp`  generalized-potential-enstrophy mode EVP

## Discussion

  At fixed positive horizontal wavenumber `k`, this factory
  creates the `F` problem
  $$
  -\frac{\partial}{\partial z}\left(\frac{f_0^2}{N^2}\frac{\partial F_j}{\partial z}\right)
  +k^2F_j=\Lambda_jF_j.
  $$
  A finite surface weight applies
  $$
  \frac{f_0^2}{N_s^2}F_j'(z_s)+\frac{f_0^2}{g}F_j(z_s)
  =\Lambda_j\frac{f_0^2}{\alpha_0}F_j(z_s),
  $$
  while a finite bottom weight applies
  $$
  \frac{f_0^2}{N_b^2}F_j'(z_b)
  =-\Lambda_j\frac{f_0^2}{\alpha_d}F_j(z_b).
  $$
  Positive infinity removes the corresponding eigenvalue-side
  endpoint term. The default bottom is inactive. When `alpha0`
  is omitted, the surface weight is
  $$
  \alpha_0=\frac{f_0^2}{b_\mathrm{eff}},\qquad
  b_\mathrm{eff}=\frac{(\int N\,dz)^2}{4\int N^2\,dz}.
  $$

  Solved bases use
  `Normalization.generalizedPotentialEnstrophy` by default, so
  $$
  \frac{1}{D}\left[\int F_iF_j\,dz
  +\frac{f_0^2}{\alpha_0}F_i(z_s)F_j(z_s)
  +\frac{f_0^2}{\alpha_d}F_i(z_b)F_j(z_b)\right]=\delta_{ij},
  $$
  with inactive terms omitted. The equivalent depth is the
  diagnostic quantity
  $$h_j=f_0^2/[g(\Lambda_j-k^2)].$$

  ```matlab
  evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
      N2=N2,zDomain=[-4000 0],k=2*pi/100e3,f0=1e-4);
  basisSet = IMSolverSpectral(nEVP=128,coordinateKind="wkb").solveEVP(evp,nModes=8);
  ```
