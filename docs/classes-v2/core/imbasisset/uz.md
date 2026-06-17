---
layout: default
title: uz
parent: IMBasisSet
grand_parent: Core
nav_order: 24
mathjax: true
---

#  uz

Evaluate solved scalar vertical derivatives.


---

## Declaration
```matlab
 values = uz(basisSet,z,options)
```
## Parameters
+ `z`  physical coordinate
+ `options.normalization`  normalization to apply

## Returns
+ `values`  derivative mode values

## Discussion

  Derivatives are scaled by the same modal factors used for
  `u`, so $$u'_j(z)=u_j^{\mathrm{raw}\prime}(z)/s_j$$.
