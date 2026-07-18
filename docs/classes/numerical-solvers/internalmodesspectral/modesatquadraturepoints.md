---
layout: default
title: modesAtQuadraturePoints
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 56
mathjax: true
---

#  modesAtQuadraturePoints

Return resolved modes on mode-adapted quadrature points.


---

## Declaration
```matlab
 [F,G,h,zq] = modesAtQuadraturePoints(self,options)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `options.nPoints`  number of quadrature points and modes requested
+ `options.omega`  fixed frequency in radians per second; specify either `omega` or `k`
+ `options.k`  fixed horizontal wavenumber; specify either `omega` or `k`

## Returns
+ `F`  horizontal-velocity mode matrix on `zq`
+ `G`  vertical-velocity mode matrix on `zq`
+ `h`  equivalent-depth row vector
+ `zq`  quadrature-point depth locations

## Discussion

This helper chooses quadrature points from the roots of a resolved
higher mode, temporarily evaluates the requested modes on those points,
and restores the original output grid before returning.
