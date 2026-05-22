---
layout: default
title: F
parent: IMBasisSet
grand_parent: Core
nav_order: 1
mathjax: true
---

#  F

Evaluate $$F_j(z)$$ modes on a physical grid.


---

## Declaration
```matlab
 F = F(basisSet,z,options)
```
## Parameters
+ `z`  physical evaluation points
+ `options.normalization`  normalization to apply

## Returns
+ `F`  evaluated `F` modes

## Discussion

  If the EVP formulation solves `F`, this returns the solved
  modes. If the EVP formulation solves `G`, this evaluates
  $$F_j=h_j\partial_zG_j$$.
