---
layout: default
title: G
parent: IMBasisSet
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  G

Evaluate $$G_j(z)$$ modes on a physical grid.


---

## Declaration
```matlab
 G = G(basisSet,z,options)
```
## Parameters
+ `z`  physical evaluation points
+ `options.normalization`  normalization to apply

## Returns
+ `G`  evaluated `G` modes

## Discussion

  If the EVP formulation solves `G`, this returns the solved
  modes. If the EVP formulation solves `F`, this evaluates
  $$G_j=-gN^{-2}\partial_zF_j$$.
