---
layout: default
title: geostrophicNormFactor
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 10
mathjax: true
---

#  geostrophicNormFactor

Return the geostrophic normalization factor.

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

This developer utility implements the analytical
`modeFamily="hydrostatic"` normalization rule. Baroclinic
modes use one shared scale factor for the coupled `F`/`G`
pair based on the raw `G` inner product. A barotropic zero
mode uses the `F` norm divided by
$$\sqrt{z_\mathrm{surface}-z_\mathrm{bottom}}$$.
