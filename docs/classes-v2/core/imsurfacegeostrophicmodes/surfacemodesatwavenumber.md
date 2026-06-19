---
layout: default
title: surfaceModesAtWavenumber
parent: IMSurfaceGeostrophicModes
grand_parent: Core
nav_order: 9
mathjax: true
---

#  surfaceModesAtWavenumber

Create surface SQG modes at fixed wavenumber.


---

## Declaration
```matlab
 problem = IMSurfaceGeostrophicModes.surfaceModesAtWavenumber(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.f0`  Coriolis parameter
+ `options.k`  horizontal wavenumbers
+ `options.metadata`  additional metadata

## Returns
+ `problem`  surface SQG problem

## Discussion
