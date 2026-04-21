---
layout: default
title: DistributePointsInRegionsWithMinimum
parent: InternalModesAdaptiveSpectral
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  DistributePointsInRegionsWithMinimum

Distribute collocation points across adaptive regions.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 nEVPPoints = DistributePointsInRegionsWithMinimum(nTotalPoints,boundaries,thesign)
```
## Parameters
+ `nTotalPoints`  total number of collocation points
+ `boundaries`  region boundaries in WKB coordinates
+ `thesign`  signs of `N2-\omega^2` in each region

## Returns
+ `nEVPPoints`  point counts assigned to each region

## Discussion

                This algorithm distributes the points/polynomials amongst the
  different regions/equations
