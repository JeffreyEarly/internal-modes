---
layout: default
title: normalizationNames
parent: IMBasisSet
grand_parent: Core
nav_order: 17
mathjax: true
---

#  normalizationNames

Return installed normalization rule names.


---

## Declaration
```matlab
 names = normalizationNames(basisSet)
```
## Returns
+ `names`  string array of installed normalization names

## Discussion

`normalizationNames` reports the rules available to
`normalizationFactors` and selectable by
`basisSet.normalization`.
