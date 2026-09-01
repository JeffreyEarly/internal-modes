---
layout: default
title: crossSpectrum
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 6
mathjax: true
---

#  crossSpectrum

Compute an internal-mode modal cross-spectrum.


---

## Declaration
```matlab
 spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB,options)
```
## Parameters
+ `coefficientsA`  first modal coefficients
+ `coefficientsB`  second modal coefficients
+ `options.variable`  optional variable name, `"F"` or `"G"`

## Returns
+ `spectrum`  modal cross-spectrum

## Discussion

  If `options.variable` is omitted, the solved formulation is
  used. The requested variable must have a known inner product.
