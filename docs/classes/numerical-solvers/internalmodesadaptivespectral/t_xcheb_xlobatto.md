---
layout: default
title: T_xCheb_xLobatto
parent: InternalModesAdaptiveSpectral
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  T_xCheb_xLobatto

Transform coupled Chebyshev coefficients onto the adaptive Lobatto grid.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 v_xLobatto = T_xCheb_xLobatto(self,v_xCheb)
```
## Parameters
+ `self`  InternalModesAdaptiveSpectral instance
+ `v_xCheb`  coefficients in the coupled Chebyshev basis

## Returns
+ `v_xLobatto`  values on the adaptive Lobatto grid

## Discussion

              transform from xCheb basis to xLobatto
