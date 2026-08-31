---
layout: default
title: rotateWithPencil
parent: IMAnalyticalGeostrophicZeroAPVModesBasis
grand_parent: Analytical bases
nav_order: 19
mathjax: true
---

#  rotateWithPencil

Apply a custom symmetric matrix-pencil rotation.


---

## Declaration
```matlab
 basisSet = rotateWithPencil(exactModes,options)
```
## Parameters
+ `options.name`  custom rotation name
+ `options.leftMatrix`  symmetric left matrix pages
+ `options.rightMatrix`  symmetric right matrix pages

## Returns
+ `basisSet`  custom-pencil rotated exact basis

## Discussion

  Input matrices are expressed in canonical endpoint coordinates.
  Two-dimensional matrices broadcast across wavenumber pages.
