---
layout: default
title: normalizedNativeModes
parent: IMBasisSet
grand_parent: Core
nav_order: 18
mathjax: true
---

#  normalizedNativeModes

Return native modes scaled by a normalization.


---

## Declaration
```matlab
 modes = normalizedNativeModes(basisSet,normalization)
```
## Parameters
+ `normalization`  normalization convention

## Returns
+ `modes`  scaled native mode columns

## Discussion

  The native coefficient columns are divided by the same factors
  returned by `normalizationFactors`.
