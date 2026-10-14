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
+ `options.rotationRate`  planetary rotation rate in radians per second; default Earth rotation
+ `options.latitude`  central latitude in degrees north; default 33
+ `options.g`  gravitational acceleration in meters per second squared
+ `options.alpha0`  optional positive surface generalized-potential-enstrophy weight
+ `options.alphaD`  positive bottom generalized-potential-enstrophy weight or `Inf`

## Returns
+ `evp`  generalized-potential-enstrophy mode EVP

## Discussion

At fixed positive horizontal wavenumber `k`, numerical solves
use the regular pair
$$F_j'=-N^2G_j/g,\qquad G_j'=F_j/h_j,$$
with public eigenvalue
$$\Lambda_j=k^2+f_0^2/(g h_j).$$
The endpoint conditions are
$$\alpha_0(G_j(z_s)-F_j(z_s))=-g\Lambda_jF_j(z_s),$$
$$\alpha_DG_j(z_b)=g\Lambda_jF_j(z_b).$$
This formulation allows finite real nonnegative N2. Its
inherited scalar coefficients apply where N2 is positive.
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
with inactive terms omitted. The equivalent depth is retained
from the solved inverse depth, avoiding subtraction of nearly
equal public eigenvalues; exact zero gives h=Inf. Public
eigenvalues use a positive Rayleigh quotient. Fully neutral
stratification requires an explicit alpha0 because its default
stratification scale is undefined.

```matlab
evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
    N2=N2,zDomain=[-4000 0],k=2*pi/100e3,rotationRate=7.2921e-5,latitude=45);
basisSet = IMSolverSpectral(nEVP=128,coordinateKind="wkb").solveEVP(evp,nModes=8);
```
