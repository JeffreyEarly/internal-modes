---
layout: default
title: hydrostaticGModes
parent: IMInternalModes
grand_parent: Core
nav_order: 8
mathjax: true
---

#  hydrostaticGModes

Create the hydrostatic `G` internal-mode EVP.


---

## Declaration
```matlab
 evp = IMInternalModes.hydrostaticGModes(options)
```
## Parameters
+ `options.N2`  buoyancy frequency squared function
+ `options.zDomain`  physical vertical domain
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.surfaceBoundary`  surface endpoint condition
+ `options.bottomBoundary`  bottom endpoint condition

## Returns
+ `evp`  hydrostatic `G` EVP

## Discussion

  The canonical scalar form is
  $$-G''=\lambda N^2G/g.$$
