---
layout: default
title: projectOntoGModesAtWavenumber
parent: InternalModes
grand_parent: Classes
nav_order: 24
mathjax: true
---

#  projectOntoGModesAtWavenumber

Project a profile onto the `G` modes at fixed horizontal wavenumber.


---

## Declaration
```matlab
 [m,G] = projectOntoGModesAtWavenumber(self,zeta,k)
```
## Parameters
+ `self`  InternalModes instance
+ `zeta`  profile sampled on `zOut`
+ `k`  horizontal wavenumber

## Returns
+ `m`  modal coefficients from the least-squares projection
+ `G`  `G`-mode matrix used in the projection

## Discussion


