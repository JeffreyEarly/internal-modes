---
layout: default
title: FindTurningPointBoundariesAtFrequency
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 20
mathjax: true
---

#  FindTurningPointBoundariesAtFrequency

Find turning points and region signs for $$N^2(z)-\omega^2$$.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [zBoundariesAndTPs,thesign,boundaryIndices] = FindTurningPointBoundariesAtFrequency(N2,z,omega)
```
## Parameters
+ `N2`  buoyancy-frequency profile on grid `z`
+ `z`  depth grid
+ `omega`  target frequency in radians per second

## Returns
+ `zBoundariesAndTPs`  boundaries consisting of the endpoints and turning points
+ `thesign`  sign of $$N^2-\omega^2$$ in each region
+ `boundaryIndices`  grid indices associated with each returned boundary

## Discussion

                    This function returns not just the turning points, but also
  the top and bottom boundary locations in z. The boundary
  indices are the index to the point just *above* the turning
  point.
