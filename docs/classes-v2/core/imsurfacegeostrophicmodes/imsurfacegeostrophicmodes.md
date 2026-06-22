---
layout: default
title: IMSurfaceGeostrophicModes
parent: IMSurfaceGeostrophicModes
grand_parent: Core
nav_order: 1
mathjax: true
---

#  IMSurfaceGeostrophicModes

Create a projected surface-geostrophic boundary-mode problem.


---

## Declaration
```matlab
 problem = IMSurfaceGeostrophicModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.k`  horizontal wavenumbers
+ `options.g0`  surface buoyancy-anomaly weight
+ `options.gd`  bottom buoyancy-anomaly weight
+ `options.surfaceAnomaly`  surface anomaly convention
+ `options.metadata`  additional metadata

## Returns
+ `problem`  surface-geostrophic boundary-mode problem

## Discussion
