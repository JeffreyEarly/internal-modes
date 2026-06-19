---
layout: default
title: bottomModesAtWavenumber
parent: IMSurfaceGeostrophicModes
grand_parent: Core
nav_order: 3
mathjax: true
---

#  bottomModesAtWavenumber

Create bottom SQG modes at fixed wavenumber.


---

## Declaration
```matlab
 problem = IMSurfaceGeostrophicModes.bottomModesAtWavenumber(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.f0`  Coriolis parameter
+ `options.k`  horizontal wavenumbers
+ `options.metadata`  additional metadata

## Returns
+ `problem`  bottom SQG problem

## Discussion
