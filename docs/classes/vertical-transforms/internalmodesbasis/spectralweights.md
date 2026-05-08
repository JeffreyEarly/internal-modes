---
layout: default
title: spectralWeights
parent: InternalModesBasis
grand_parent: Classes
nav_order: 29
mathjax: true
---

#  spectralWeights

Return Parseval weights for modal spectra.


---

## Declaration
```matlab
 weights = spectralWeights(self,component,options)
```
## Parameters
+ `self`  InternalModesBasis instance
+ `component`  `"F"` or `"G"`
+ `options.nModes`  number of weights requested

## Returns
+ `weights`  column vector of modal spectrum weights

## Discussion

  For canonical geostrophic F modes this returns
  $$\gamma_0=D$$ and $$\gamma_j=h_g^j$$. For canonical G modes
  the current normalization gives a modal spectrum weighted by
  $$g$$.
