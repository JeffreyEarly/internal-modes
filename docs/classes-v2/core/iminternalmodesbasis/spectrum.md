---
layout: default
title: spectrum
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 14
mathjax: true
---

#  spectrum

Compute an internal-mode modal spectrum.


---

## Declaration
```matlab
 spectrum = spectrum(basisSet,coefficients,options)
```
## Parameters
+ `coefficients`  modal coefficients
+ `options.variable`  optional variable name, `"F"` or `"G"`

## Returns
+ `spectrum`  modal spectrum

## Discussion

  If `options.variable` is omitted, the solved formulation is
  used. The requested variable must have a known inner product.
