---
layout: default
title: ProjectOntoGModesAtWavenumber
parent: InternalModes
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  ProjectOntoGModesAtWavenumber

Project a profile onto the `G` modes at fixed horizontal wavenumber.


---

## Declaration
```matlab
 [m,G] = ProjectOntoGModesAtWavenumber(self,zeta,k)
```
## Parameters
+ `self`  InternalModes instance
+ `zeta`  profile sampled on `zOut`
+ `k`  horizontal wavenumber

## Returns
+ `m`  modal coefficients from the least-squares projection
+ `G`  `G`-mode matrix used in the projection

## Discussion


