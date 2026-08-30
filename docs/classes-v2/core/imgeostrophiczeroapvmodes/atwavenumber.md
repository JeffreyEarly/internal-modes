---
layout: default
title: atWavenumber
parent: IMGeostrophicZeroAPVModes
grand_parent: Core
nav_order: 3
mathjax: true
---

#  atWavenumber

Create geostrophic zero-APV modes at fixed wavenumber.


---

## Declaration
```matlab
 problem = IMGeostrophicZeroAPVModes.atWavenumber(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.k`  horizontal wavenumbers
+ `options.endpoints`  requested endpoint coordinates
+ `options.surfaceBoundary`  surface endpoint convention
+ `options.metadata`  additional metadata

## Returns
+ `problem`  geostrophic zero-APV problem

## Discussion
