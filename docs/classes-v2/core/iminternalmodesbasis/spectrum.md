---
layout: default
title: spectrum
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 23
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
  This is a signed spectrum: entries associated with negative
  Pontryagin directions can be negative. Use `majorantNorm` for
  a positive total magnitude; a generally additive per-mode
  majorant spectrum does not exist because the majorant Gram
  matrix need not be diagonal in this basis.
