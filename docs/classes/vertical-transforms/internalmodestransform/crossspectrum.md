---
layout: default
title: crossSpectrum
parent: InternalModesTransform
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  crossSpectrum

Compute a vertical modal cross-spectrum.


---

## Declaration
```matlab
 S = crossSpectrum(self,coefficientsA,coefficientsB,options)
```
## Parameters
+ `self`  InternalModesTransform instance
+ `coefficientsA`  first modal coefficient array
+ `coefficientsB`  second modal coefficient array
+ `options.component`  `"F"` or `"G"`
+ `options.horizontalMultiplicity`  multiplicity factor, often 1 or 2
+ `options.requireCanonical`  true to reject noncanonical component spectra

## Returns
+ `S`  modal cross-spectrum

## Discussion

  For component spectral weights $$s_j$$, coefficients $$a_j$$
  and $$b_j$$, and multiplicity $$m$$, this method returns

  $$
  S_j = m\,s_j\,\Re\{a_j b_j^*\}.
  $$
