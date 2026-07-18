---
layout: default
title: u
parent: IMBasisSet
grand_parent: Core
nav_order: 26
mathjax: true
---

#  u

Evaluate solved scalar modes.


---

## Declaration
```matlab
 values = u(basisSet,z,options)
```
## Parameters
+ `z`  physical coordinate
+ `options.normalization`  normalization to apply

## Returns
+ `values`  scalar mode values

## Discussion

Returned columns are normalized as
$$u_j(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
where $$s_j$$ comes from `normalizationFactors`.
