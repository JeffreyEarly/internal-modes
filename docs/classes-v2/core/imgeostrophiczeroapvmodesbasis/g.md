---
layout: default
title: G
parent: IMGeostrophicZeroAPVModesBasis
grand_parent: Core
nav_order: 2
mathjax: true
---

#  G

Evaluate diagnostic displacement structures $$G(z)$$.


---

## Declaration
```matlab
 values = G(basisSet,z)
```
## Parameters
+ `z`  physical coordinate

## Returns
+ `values`  page-shaped `G` values

## Discussion

The result has dimensions `nZ x nEndpoints x nK` and uses
$$G=-gN^{-2}\partial_zF$$.
