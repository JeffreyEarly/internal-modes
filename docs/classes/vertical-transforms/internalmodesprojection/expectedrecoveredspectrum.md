---
layout: default
title: expectedRecoveredSpectrum
parent: InternalModesProjection
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  expectedRecoveredSpectrum

Apply the spectral window to a true candidate-mode spectrum.


---

## Declaration
```matlab
 expectedSpectrum = expectedRecoveredSpectrum(self,trueSpectrum)
```
## Parameters
+ `self`  InternalModesProjection instance
+ `trueSpectrum`  candidate-mode spectrum

## Returns
+ `expectedSpectrum`  expected retained recovered spectrum

## Discussion

For uncorrelated modal coefficients, `spectralWindow` maps
the candidate spectrum into the expected retained recovered
spectrum.
