---
layout: default
title: addNormalization
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 5
mathjax: true
---

#  addNormalization

Add a named normalization rule.


---

## Declaration
```matlab
 basisSet = addNormalization(basisSet,name,rule)
```
## Parameters
+ `name`  normalization rule name
+ `rule`  function handle with signature `scale = rule(basisSet,iMode)`

## Returns
+ `basisSet`  basis set with the rule installed

## Discussion
