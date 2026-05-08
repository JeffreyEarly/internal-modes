---
layout: default
title: crossTransformTo
parent: InternalModesTransform
grand_parent: Classes
nav_order: 10
mathjax: true
---

#  crossTransformTo

Map coefficients from a target transform into this transform.


---

## Declaration
```matlab
 matrix = crossTransformTo(self,targetTransform,options)
```
## Parameters
+ `self`  InternalModesTransform instance
+ `targetTransform`  transform whose coefficients should be evaluated in this transform
+ `options.component`  `"F"` or `"G"`
+ `options.allowNoncanonical`  true to allow numerical F matrices

## Returns
+ `matrix`  coefficient-to-coefficient transform matrix

## Discussion

  The matrix is

  $$
  C = P_{\mathrm{self}}\Phi_{\mathrm{target}},
  $$

  where `P_self` is this transform's forward matrix and
  `Phi_target` is the target transform's inverse matrix for the
  selected component.
