---
layout: default
title: spectrum
parent: InternalModesProjection
grand_parent: Classes
nav_order: 29
mathjax: true
---

#  spectrum

Compute a retained-mode observation spectrum.


---

## Declaration
```matlab
 S = spectrum(self,coefficients,options)
```
## Parameters
+ `self`  InternalModesProjection instance
+ `coefficients`  retained modal coefficients
+ `options.horizontalMultiplicity`  multiplicity factor, often 1 or 2
+ `options.requireCanonical`  true to reject noncanonical spectra

## Returns
+ `S`  retained-mode spectrum

## Discussion

  This method uses the retained rows of the candidate spectral
  weights. For canonical G modes this is the potential-energy
  spectrum implied by the basis normalization.
