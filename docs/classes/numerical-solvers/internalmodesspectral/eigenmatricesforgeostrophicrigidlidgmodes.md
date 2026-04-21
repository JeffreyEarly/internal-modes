---
layout: default
title: EigenmatricesForGeostrophicRigidLidGModes
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  EigenmatricesForGeostrophicRigidLidGModes

Assemble the rigid-lid geostrophic EVP with displaced boundaries.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [A,B] = EigenmatricesForGeostrophicRigidLidGModes(self,eta0,etad)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `eta0`  imposed surface displacement coefficient
+ `etad`  imposed bottom displacement coefficient

## Returns
+ `A`  left generalized-eigenproblem matrix
+ `B`  right generalized-eigenproblem matrix

## Discussion


