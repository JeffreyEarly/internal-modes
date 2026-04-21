---
layout: default
title: GaussQuadraturePointsForModesAtFrequency
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 26
mathjax: true
---

#  GaussQuadraturePointsForModesAtFrequency

Return quadrature points tailored to fixed-$$\omega$$ modes.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 z_g = GaussQuadraturePointsForModesAtFrequency(self,nPoints,omega)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `nPoints`  number of quadrature points requested
+ `omega`  target frequency in radians per second

## Returns
+ `z_g`  depth locations of the quadrature points

## Discussion


