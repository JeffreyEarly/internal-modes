---
layout: default
title: rawU
parent: IMBasisSet
grand_parent: Core
nav_order: 19
mathjax: true
---

#  rawU

Evaluate raw solved scalar modes.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 values = rawU(basisSet,z)
```
## Parameters
+ `z`  physical coordinate

## Returns
+ `values`  unnormalized scalar mode values

## Discussion

  This developer utility evaluates the unnormalized native
  mode columns on the physical grid `z`. The public `u` method
  applies the active normalization rule after this step:
  $$u_j(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
  where $$s_j$$ comes from `normalizationFactors`.
