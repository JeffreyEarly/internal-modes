---
layout: default
title: evaluate
parent: IMBasisSet
grand_parent: Core
nav_order: 4
mathjax: true
---

#  evaluate

the scalar variable.


---

## Declaration
```matlab
 values = evaluate(basisSet,variable,z,options)
```
## Parameters
+ `variable`  `"u"`
+ `z`  physical coordinate
+ `options.normalization`  normalization to apply

## Returns
+ `values`  scalar mode values

## Discussion

  `evaluate("u",z)` is equivalent to `u(z)` and accepts the
  same `normalization` override.
