---
layout: default
title: ModesAtFrequency
parent: InternalModesBase
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  ModesAtFrequency

Return vertical modes for a fixed frequency.


---

## Declaration
```matlab
 [F,G,h,k] = ModesAtFrequency(self,omega)
```
## Parameters
+ `self`  concrete InternalModesBase subclass instance
+ `omega`  frequency in radians per second

## Returns
+ `F`  horizontal-velocity mode matrix on `zOut`
+ `G`  vertical-velocity mode matrix on `zOut`
+ `h`  equivalent-depth row vector
+ `k`  horizontal wavenumber row vector implied by `h` and `omega`

## Discussion


