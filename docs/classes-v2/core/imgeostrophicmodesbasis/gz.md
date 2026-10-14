---
layout: default
title: Gz
parent: IMGeostrophicModesBasis
grand_parent: Core
nav_order: 1
mathjax: true
---

#  Gz

Evaluate derivatives of the stored companion `G` modes.


---

## Declaration
```matlab
 Gz = Gz(basisSet,z,options)
```
## Parameters
+ `z`  physical coordinate
+ `options.normalization`  normalization rule name or enum value

## Returns
+ `Gz`  evaluated `G` derivatives

## Discussion

This differentiates the native representation independently of
the first-order differential equations.
