---
layout: default
title: geostrophicNormFactor
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 9
mathjax: true
---

#  geostrophicNormFactor

Return the hydrostatic geostrophic normalization factor.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 factor = geostrophicNormFactor(basisSet,iMode)
```
## Parameters
+ `iMode`  retained mode index

## Returns
+ `factor`  raw geostrophic scale factor

## Discussion

This developer utility implements the normalization rule
installed for `modeFamily="hydrostatic"`. It chooses one
shared raw scale factor for each coupled `F`/`G` mode.
Baroclinic modes use
$$s_j^2=\langle G_j,G_j\rangle_G,$$
so normalized modes satisfy
$$\langle G_j,G_j\rangle_G=1,\qquad
\langle F_j,F_j\rangle_F=h_j.$$
A barotropic zero mode uses the `F` norm divided by
$$\sqrt{z_\mathrm{surface}-z_\mathrm{bottom}}$$ and is a
separate null-mode convention.
