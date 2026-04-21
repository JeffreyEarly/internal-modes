---
layout: default
title: ModesAtWavenumber
parent: InternalModesBase
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  ModesAtWavenumber

Return vertical modes for a fixed horizontal wavenumber.


---

## Declaration
```matlab
 [F,G,h,omega] = ModesAtWavenumber(self,k)
```
## Parameters
+ `self`  concrete InternalModesBase subclass instance
+ `k`  horizontal wavenumber

## Returns
+ `F`  horizontal-velocity mode matrix on `zOut`
+ `G`  vertical-velocity mode matrix on `zOut`
+ `h`  equivalent-depth row vector
+ `omega`  frequency row vector implied by `h` and `k`

## Discussion


