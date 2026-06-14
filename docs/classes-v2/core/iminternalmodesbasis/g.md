---
layout: default
title: G
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 2
mathjax: true
---

#  G

Evaluate `G` modes.


---

## Declaration
```matlab
 G = G(basisSet,z,options)
```
## Parameters
+ `z`  physical coordinate
+ `options.normalization`  normalization to apply

## Returns
+ `G`  evaluated `G` modes

## Discussion

  If the EVP formulation is `G`, this evaluates the solved
  canonical variable. If the formulation is `F`, `G` is recovered
  from $$G=-gN^{-2}F_z$$.
