---
layout: default
title: waveModesAtFrequency
parent: IMInternalModes
grand_parent: Core
nav_order: 7
mathjax: true
---

#  waveModesAtFrequency

Create the fixed-frequency wave-mode EVP.


---

## Declaration
```matlab
 evp = IMInternalModes.waveModesAtFrequency(options)
```
## Parameters
+ `options.omega`  wave frequency
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface endpoint condition
+ `options.bottomBoundary`  bottom endpoint condition

## Returns
+ `evp`  fixed-frequency `G` EVP

## Discussion

  The canonical scalar form is
  $$-G''=\lambda(N^2-\omega^2)G/g.$$
