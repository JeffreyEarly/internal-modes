---
layout: default
title: projectVertical
parent: InternalModesTransform
grand_parent: Classes
nav_order: 31
mathjax: true
---

#  projectVertical

Project vertical samples onto modal coefficients.


---

## Declaration
```matlab
 coefficients = projectVertical(self,field,options)
```
## Parameters
+ `self`  InternalModesTransform instance
+ `field`  vertical samples with rows matching `z`
+ `options.component`  `"F"` or `"G"`
+ `options.allowNoncanonical`  true to allow numerical wave-F projection

## Returns
+ `coefficients`  modal coefficients

## Discussion

  If `field` has multiple columns, each column is projected
  independently with the same vertical operator.
