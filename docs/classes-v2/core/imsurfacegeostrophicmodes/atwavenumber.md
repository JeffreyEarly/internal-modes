---
layout: default
title: atWavenumber
parent: IMSurfaceGeostrophicModes
grand_parent: Core
nav_order: 3
mathjax: true
---

#  atWavenumber

Create projected surface-geostrophic modes at fixed wavenumber.


---

## Declaration
```matlab
 problem = IMSurfaceGeostrophicModes.atWavenumber(options)
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
+ `problem`  projected SQG problem

## Discussion
