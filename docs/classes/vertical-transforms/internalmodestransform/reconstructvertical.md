---
layout: default
title: reconstructVertical
parent: InternalModesTransform
grand_parent: Classes
nav_order: 33
mathjax: true
---

#  reconstructVertical

Reconstruct vertical samples from modal coefficients.


---

## Declaration
```matlab
 field = reconstructVertical(self,coefficients,options)
```
## Parameters
+ `self`  InternalModesTransform instance
+ `coefficients`  modal coefficients
+ `options.component`  `"F"` or `"G"`

## Returns
+ `field`  reconstructed vertical samples on `z`

## Discussion
