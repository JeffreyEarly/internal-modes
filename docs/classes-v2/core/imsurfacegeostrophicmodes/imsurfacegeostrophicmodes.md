---
layout: default
title: IMSurfaceGeostrophicModes
parent: IMSurfaceGeostrophicModes
grand_parent: Core
nav_order: 1
mathjax: true
---

#  IMSurfaceGeostrophicModes

Create a surface-geostrophic boundary-mode problem.


---

## Declaration
```matlab
 problem = IMSurfaceGeostrophicModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.f0`  Coriolis parameter
+ `options.k`  horizontal wavenumbers
+ `options.boundary`  active boundary, `"surface"` or `"bottom"`
+ `options.metadata`  additional metadata

## Returns
+ `problem`  surface-geostrophic boundary-mode problem

## Discussion
