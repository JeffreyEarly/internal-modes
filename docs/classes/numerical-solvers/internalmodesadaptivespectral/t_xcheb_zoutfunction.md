---
layout: default
title: T_xCheb_zOutFunction
parent: InternalModesAdaptiveSpectral
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  T_xCheb_zOutFunction

Transform coupled Chebyshev coefficients onto the public output grid.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 v_zOut = T_xCheb_zOutFunction(self,v_xCheb)
```
## Parameters
+ `self`  InternalModesAdaptiveSpectral instance
+ `v_xCheb`  coefficients in the coupled Chebyshev basis

## Returns
+ `v_zOut`  values on `zOut`

## Discussion

              transform from xCheb basis to zOut
