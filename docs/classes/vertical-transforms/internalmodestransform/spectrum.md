---
layout: default
title: spectrum
parent: InternalModesTransform
grand_parent: Classes
nav_order: 42
mathjax: true
---

#  spectrum

Compute a vertical modal auto-spectrum.


---

## Declaration
```matlab
 S = spectrum(self,coefficients,options)
```
## Parameters
+ `self`  InternalModesTransform instance
+ `coefficients`  vertical modal coefficients
+ `options.component`  `"F"` or `"G"`
+ `options.horizontalMultiplicity`  multiplicity factor, often 1 or 2
+ `options.requireCanonical`  true to reject noncanonical component spectra

## Returns
+ `S`  modal auto-spectrum

## Discussion

This is equivalent to `crossSpectrum(coefficients,coefficients)`.
