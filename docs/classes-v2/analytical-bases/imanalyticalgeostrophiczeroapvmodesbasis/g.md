---
layout: default
title: G
parent: IMAnalyticalGeostrophicZeroAPVModesBasis
grand_parent: Analytical bases
nav_order: 2
mathjax: true
---

#  G

Evaluate exact diagnostic displacement structures $$G(z)$$.


---

## Declaration
```matlab
 values = G(basisSet,z,options)
```
## Parameters
+ `z`  physical coordinate
+ `options.pages`  source wavenumber pages, preserving order and repeats

## Returns
+ `values`  page-shaped exact `G` values

## Discussion

  The evaluator uses the closed-form derivative relation
  $$G=-gN^{-2}\partial_zF$$ without numerical differentiation.
