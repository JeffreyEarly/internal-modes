---
layout: default
title: geostrophicAPVModes
parent: IMInternalModes
grand_parent: Core
nav_order: 9
mathjax: true
---

#  geostrophicAPVModes

Create signed generalized-energy geostrophic APV modes.


---

## Declaration
```matlab
 evp = IMInternalModes.geostrophicAPVModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function in radians squared per second squared
+ `options.zDomain`  physical vertical domain in meters
+ `options.g`  gravitational acceleration in meters per second squared
+ `options.g0`  signed surface acceleration in meters per second squared
+ `options.gd`  signed bottom acceleration in meters per second squared
+ `options.surfaceBoundary`  `"freeSurface"` or `"rigidLid"`

## Returns
+ `evp`  generalized-energy geostrophic APV EVP

## Discussion

This factory creates the hydrostatic `F` problem
$$-\frac{\partial}{\partial z}\left(\frac{1}{N^2}\frac{\partial F_j}{\partial z}\right)=\frac{F_j}{g h_j}$$
with diagnostic relation
$$G_j=-\frac{g}{N^2}\frac{\partial F_j}{\partial z}.$$
The endpoint accelerations `g0` and `gd` are required and
may be signed finite values, zero, or positive infinity.
Zero selects Dirichlet; positive infinity selects the
corresponding reciprocal-zero endpoint limit.

The surface convention maps to canonical ordinary boundary
conditions as follows:

| Convention | `g0` | Canonical `(a,b)` | Physical condition |
| --- | --- | --- | --- |
| free surface | finite nonzero | $$(-(1/g+1/g_0),1)$$ | $$G_s=(1+g/g_0)F_s$$ |
| free surface | zero | Dirichlet | $$F_s=0$$ |
| free surface | `Inf` | $$(-1/g,1)$$ | $$G_s=F_s$$ |
| rigid lid | finite nonzero | $$(-1/g_0,1)$$ | $$G_s=(g/g_0)F_s$$ |
| rigid lid | zero | Dirichlet | $$F_s=0$$ |
| rigid lid | `Inf` | Neumann | $$G_s=0$$ |

The bottom condition is:

| `gd` | Canonical `(a,b)` | Physical condition |
| --- | --- | --- |
| finite nonzero | $$(1/g_d,1)$$ | $$G_b=-(g/g_d)F_b$$ |
| zero | Dirichlet | $$F_b=0$$ |
| `Inf` | Neumann | $$G_b=0$$ |

The signed generalized-energy inner products are
$$\langle F_i,F_j\rangle_F=\int_{z_b}^{z_s}F_iF_j\,dz$$
and
$$\langle G_i,G_j\rangle_G=\int_{z_b}^{z_s}\frac{N^2}{g}G_iG_j\,dz+c_sG_i(z_s)G_j(z_s)+\frac{g_d}{g}G_i(z_b)G_j(z_b).$$
For finite active endpoints, $$c_s=g_0/(g+g_0)$$
under the free-surface convention and $$c_s=g_0/g$$
under the rigid-lid convention. Dirichlet, Neumann,
positive-infinite, and the free-surface `g0=-g` reductions
omit any coefficient made singular by the reduced endpoint
constraint. The free-surface `g0=Inf` limit instead has
$$c_s=1$$.

Solved bases use `Normalization.depth` by default, so
$$D^{-1}\int_{z_b}^{z_s}F_j^2\,dz=1$$ without endpoint
terms. The same positive factor scales `F` and `G`, including
negative-eigendepth modes. Signed `h` values are retained;
an exact zero eigenvalue is represented by `h=Inf`.
`g0`, `gd`, and the surface convention are copied into
`evp.parameters` and basis metadata.
`surfaceBoundary` is the endpoint-convention contract used
to identify the matching zero-APV family. `N2(z)` must
return values with the shape of `z`; evaluated basis methods
return one column per retained mode.

```matlab
evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-4000 0],g0=-0.02,gd=Inf);
basisSet = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=4);
F = basisSet.F(z);
```

```matlab
evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-4000 0],g0=Inf,gd=Inf,surfaceBoundary="rigidLid");
```
