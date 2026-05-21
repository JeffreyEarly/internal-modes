---
layout: default
title: evaluate
parent: IMBasisSet
grand_parent: Classes
nav_order: 10
mathjax: true
---

#  evaluate

`F` or `G` on a physical grid.


---

## Declaration
```matlab
 values = evaluate(basisSet,variable,z,options)
```
## Parameters
+ `variable`  `"F"` or `"G"`
+ `z`  physical evaluation points
+ `options.normalization`  normalization to apply

## Returns
+ `values`  evaluated variable matrix

## Discussion

  This method is a validated dispatcher around the canonical
  `F(z)` and `G(z)` methods.
