---
layout: default
title: F
parent: IMAnalyticalGeostrophicZeroAPVModesBasis
grand_parent: Analytical bases
nav_order: 1
mathjax: true
---

#  F

Evaluate exact streamfunction structures $$F(z)$$.


---

## Declaration
```matlab
 values = F(basisSet,z,options)
```
## Parameters
+ `z`  physical coordinate
+ `options.pages`  source wavenumber pages, preserving order and repeats

## Returns
+ `values`  page-shaped exact `F` values

## Discussion

  The result has dimensions `nZ x nEndpoints x nK`.
