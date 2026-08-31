---
layout: default
title: geostrophicZeroAPVModesAtWavenumber
parent: IMConstantStratificationSolution
grand_parent: Analytical bases
nav_order: 4
mathjax: true
---

#  geostrophicZeroAPVModesAtWavenumber

Create exact canonical geostrophic zero-APV modes.


---

## Declaration
```matlab
 exactModes = geostrophicZeroAPVModesAtWavenumber(solution,k,options)
```
## Parameters
+ `k`  positive horizontal wavenumbers
+ `options.endpoints`  requested surface and bottom coordinates
+ `options.surfaceBoundary`  `"freeSurface"` or `"rigidLid"`
+ `options.metadata`  additional metadata

## Returns
+ `exactModes`  exact boundary-normalized basis

## Discussion

  The exact columns solve

  $$
  \frac{f_0^2}{N_0^2}F_{zz}-k^2F=0,
  \qquad G=-\frac{g}{N_0^2}F_z,
  $$

  and have unit response at each requested endpoint and zero
  response at the other endpoint.
