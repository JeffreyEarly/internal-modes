---
layout: default
title: WKBDecaySolution
parent: InternalModesAdaptiveSpectral
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  WKBDecaySolution

Evaluate the Airy-based WKB decay envelope used for pruning regions.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 decay = WKBDecaySolution(xi,L_osc,N2Omega2)
```
## Parameters
+ `xi`  positive distance from a turning point in WKB coordinates
+ `L_osc`  total oscillatory extent used in the asymptotic scaling
+ `N2Omega2`  values of `N2-\omega^2`

## Returns
+ `decay`  WKB decay-envelope estimate

## Discussion

                The decay part of the lowest F-mode. xi should be > 0
