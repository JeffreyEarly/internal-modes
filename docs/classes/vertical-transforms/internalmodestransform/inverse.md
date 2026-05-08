---
layout: default
title: inverse
parent: InternalModesTransform
grand_parent: Classes
nav_order: 22
mathjax: true
---

#  inverse

Return an inverse vertical reconstruction matrix.


---

## Declaration
```matlab
 matrix = inverse(self,options)
```
## Parameters
+ `self`  InternalModesTransform instance
+ `options.component`  `"F"` or `"G"`

## Returns
+ `matrix`  inverse reconstruction matrix

## Discussion

  The returned matrix maps vertical modal coefficients to field
  samples on `z`.
