---
layout: default
title: ModeFunctionsForOmegaAndC
parent: InternalModesExponentialStratification
grand_parent: Classes
nav_order: 11
mathjax: true
---

#  ModeFunctionsForOmegaAndC

Select the exact or approximate analytical mode functions.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [F,G] = ModeFunctionsForOmegaAndC(self,omega,c)
```
## Parameters
+ `self`  InternalModesExponentialStratification instance
+ `omega`  frequency row vector
+ `c`  phase-speed row vector

## Returns
+ `F`  function handle for the analytical `F(z,\omega,c)`
+ `G`  function handle for the analytical `G(z,\omega,c)`

## Discussion


