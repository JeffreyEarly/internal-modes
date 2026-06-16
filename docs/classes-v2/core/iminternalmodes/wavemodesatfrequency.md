---
layout: default
title: waveModesAtFrequency
parent: IMInternalModes
grand_parent: Core
nav_order: 16
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
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.omega`  wave frequency
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface endpoint condition
+ `options.bottomBoundary`  bottom endpoint condition

## Returns
+ `evp`  fixed-frequency `G` EVP

## Discussion

  The canonical scalar form is
  $$-\frac{\partial^2 G}{\partial z^2}(z)
  =\lambda\frac{N^2(z)-\omega^2}{g}G(z).$$
  Solved fixed-frequency basis sets use the generic `unity`
  normalization by default. A fixed-frequency diagnostic `F`
  inner-product normalization is deferred until the wave
  diagnostic inner-product catalog is derived. This factory
  adds `parameters.omega` and sets `parameters.formulation`,
  `parameters.f0`, and `parameters.g`.
