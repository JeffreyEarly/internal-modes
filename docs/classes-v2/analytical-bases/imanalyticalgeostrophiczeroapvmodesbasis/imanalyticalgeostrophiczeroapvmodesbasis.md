---
layout: default
title: IMAnalyticalGeostrophicZeroAPVModesBasis
parent: IMAnalyticalGeostrophicZeroAPVModesBasis
grand_parent: Analytical bases
nav_order: 3
mathjax: true
---

#  IMAnalyticalGeostrophicZeroAPVModesBasis

Create an exact canonical boundary-normalized basis.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 basisSet = IMAnalyticalGeostrophicZeroAPVModesBasis(options)
```
## Parameters
+ `options.solution`  analytical solution family
+ `options.problem`  canonical zero-APV problem
+ `options.FFunction`  exact canonical `F` evaluator
+ `options.GFunction`  exact canonical `G` evaluator
+ `options.metadata`  creation metadata

## Returns
+ `basisSet`  exact canonical basis

## Discussion

  Concrete analytical solution families supply exact `F` and
  `G` evaluators. Users normally construct this basis through
  `geostrophicZeroAPVModesAtWavenumber` on one of those families.
