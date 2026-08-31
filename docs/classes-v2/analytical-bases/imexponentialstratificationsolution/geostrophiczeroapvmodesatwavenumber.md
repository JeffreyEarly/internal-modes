---
layout: default
title: geostrophicZeroAPVModesAtWavenumber
parent: IMExponentialStratificationSolution
grand_parent: Analytical bases
nav_order: 5
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

  The exact modified-Bessel columns solve the zero-APV equation
  and are rescaled to unit endpoint responses under either the
  free-surface or rigid-lid convention.
