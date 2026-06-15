---
layout: default
title: F
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 1
mathjax: true
---

#  F

Evaluate `F` modes.


---

## Declaration
```matlab
 F = F(basisSet,z,options)
```
## Parameters
+ `z`  physical coordinate
+ `options.normalization`  normalization rule name or enum value

## Returns
+ `F`  evaluated `F` modes

## Discussion

  If the EVP formulation is `F`, this evaluates the solved
  canonical variable. If the formulation is `G`, `F` is recovered
  from $$F=hG_z$$.
