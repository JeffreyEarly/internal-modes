---
layout: default
title: crossSpectrum
parent: InternalModesProjection
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  crossSpectrum

Compute a retained-mode observation cross-spectrum.


---

## Declaration
```matlab
 S = crossSpectrum(self,coefficientsA,coefficientsB,options)
```
## Parameters
+ `self`  InternalModesProjection instance
+ `coefficientsA`  first retained coefficient array
+ `coefficientsB`  second retained coefficient array
+ `options.horizontalMultiplicity`  multiplicity factor, often 1 or 2
+ `options.requireCanonical`  true to reject noncanonical spectra

## Returns
+ `S`  retained-mode cross-spectrum

## Discussion
