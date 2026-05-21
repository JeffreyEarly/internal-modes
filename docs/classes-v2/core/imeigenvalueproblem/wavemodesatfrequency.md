---
layout: default
title: waveModesAtFrequency
parent: IMEigenvalueProblem
grand_parent: Classes
nav_order: 23
mathjax: true
---

#  waveModesAtFrequency

Create the wave-mode `G` EVP at fixed frequency.


---

## Declaration
```matlab
 evp = IMEigenvalueProblem.waveModesAtFrequency(options)
```
## Parameters
+ `options.omega`  fixed frequency in radians per second
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  location-free surface boundary law
+ `options.bottomBoundary`  location-free bottom boundary law

## Returns
+ `evp`  fixed-frequency wave-mode `G` EVP

## Discussion

  The physical-coordinate strong form is
  $$G_{zz}=\lambda(\omega^2-N^2)G/g$$.
