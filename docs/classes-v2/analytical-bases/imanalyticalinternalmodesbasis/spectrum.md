---
layout: default
title: spectrum
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 27
mathjax: true
---

#  spectrum

Compute a modal spectrum.


---

## Declaration
```matlab
 spectrum = spectrum(basisSet,coefficients,options)
```
## Parameters
+ `coefficients`  modal coefficients
+ `options.variable`  `"F"` or `"G"`

## Returns
+ `spectrum`  modal spectrum

## Discussion

  This is a signed spectrum and can contain negative entries in
  negative Pontryagin directions. Use `majorantNorm` for a
  positive total magnitude.
