---
layout: default
title: ApplyBoundaryConditions
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  ApplyBoundaryConditions

Apply the active surface and bottom conditions to an EVP pair.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [A,B] = ApplyBoundaryConditions(self,A,B)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `A`  left generalized-eigenproblem matrix
+ `B`  right generalized-eigenproblem matrix

## Returns
+ `A`  boundary-conditioned left matrix
+ `B`  boundary-conditioned right matrix

## Discussion


