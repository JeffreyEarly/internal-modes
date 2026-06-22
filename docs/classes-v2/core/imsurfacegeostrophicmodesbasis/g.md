---
layout: default
title: G
parent: IMSurfaceGeostrophicModesBasis
grand_parent: Core
nav_order: 2
mathjax: true
---

#  G

Evaluate diagnostic SQG displacement modes.


---

## Declaration
```matlab
 values = G(basisSet,z)
```
## Parameters
+ `z`  physical coordinate

## Returns
+ `values`  diagnostic `G` values

## Discussion

  `G(z)` is recovered from the projected `F` modes using
  $$G=-gN^{-2}\partial F/\partial z$$.
