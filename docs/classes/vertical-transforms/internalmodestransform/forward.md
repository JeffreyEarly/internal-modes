---
layout: default
title: forward
parent: InternalModesTransform
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  forward

Return a forward vertical projection matrix.


---

## Declaration
```matlab
 matrix = forward(self,options)
```
## Parameters
+ `self`  InternalModesTransform instance
+ `options.component`  `"F"` or `"G"`
+ `options.allowNoncanonical`  true to allow numerical wave-F forward matrices

## Returns
+ `matrix`  forward projection matrix

## Discussion

  The returned matrix maps samples on `z` to vertical modal
  coefficients. Noncanonical wave-F projections are returned
  only when `allowNoncanonical=true`.
