---
layout: default
title: ModesAtFrequencyAiry
parent: InternalModesWKB
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  ModesAtFrequencyAiry

Return the turning-point-aware Airy approximation for fixed $$\omega$$.


---

## Declaration
```matlab
 [F,G,h] = ModesAtFrequencyAiry(self,omega)
```
## Parameters
+ `self`  InternalModesWKB instance
+ `omega`  frequency in radians per second

## Returns
+ `F`  horizontal-velocity mode matrix on `zOut`
+ `G`  vertical-velocity mode matrix on `zOut`
+ `h`  equivalent-depth row vector

## Discussion


