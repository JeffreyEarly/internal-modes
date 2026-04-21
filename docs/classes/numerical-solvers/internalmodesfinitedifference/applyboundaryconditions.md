---
layout: default
title: ApplyBoundaryConditions
parent: InternalModesFiniteDifference
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  ApplyBoundaryConditions

Apply the current lower and upper boundary conditions to an EVP pair.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [A,B] = ApplyBoundaryConditions(self,A,B)
```
## Parameters
+ `self`  InternalModesFiniteDifference instance
+ `A`  left generalized-eigenproblem matrix
+ `B`  right generalized-eigenproblem matrix

## Returns
+ `A`  boundary-conditioned left matrix
+ `B`  boundary-conditioned right matrix

## Discussion


