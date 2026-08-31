---
layout: default
title: T_xLobatto_xCheb
parent: InternalModesAdaptiveSpectral
grand_parent: Classes
nav_order: 13
mathjax: true
---

#  T_xLobatto_xCheb

Transform adaptive-grid values into the coupled Chebyshev basis.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 v_xCheb = T_xLobatto_xCheb(self,v_xLobatto)
```
## Parameters
+ `self`  InternalModesAdaptiveSpectral instance
+ `v_xLobatto`  values on the adaptive Lobatto grid

## Returns
+ `v_xCheb`  coefficients in the coupled Chebyshev basis

## Discussion

              transform from xLobatto basis to xCheb basis
