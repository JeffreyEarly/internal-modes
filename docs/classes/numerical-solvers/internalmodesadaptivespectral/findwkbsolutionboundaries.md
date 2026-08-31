---
layout: default
title: FindWKBSolutionBoundaries
parent: InternalModesAdaptiveSpectral
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  FindWKBSolutionBoundaries

Estimate adaptive-region boundaries from turning points and WKB decay.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [zBoundariesAndTPs,thesign] = FindWKBSolutionBoundaries(N2,z,omega,requiredDecay)
```
## Parameters
+ `N2`  buoyancy-frequency profile on grid `z`
+ `z`  depth grid
+ `omega`  target frequency in radians per second
+ `requiredDecay`  minimum retained evanescent amplitude ratio

## Returns
+ `zBoundariesAndTPs`  region boundaries including turning points
+ `thesign`  sign of $$N^2-\omega^2$$ in each region

## Discussion

                    The requiredDecay is <=1 and find the point at which the
  solution has decayed below that value.
