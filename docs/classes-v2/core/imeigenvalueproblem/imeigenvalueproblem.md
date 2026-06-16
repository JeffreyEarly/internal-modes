---
layout: default
title: IMEigenvalueProblem
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 1
mathjax: true
---

#  IMEigenvalueProblem

Create a canonical scalar EVP.


---

## Declaration
```matlab
 evp = IMEigenvalueProblem(options)
```
## Parameters
+ `options.name`  short EVP name
+ `options.zDomain`  physical vertical domain
+ `options.p`  derivative-flux coefficient
+ `options.q`  left-side value coefficient
+ `options.r`  eigenvalue-side metric coefficient
+ `options.surfaceBoundary`  surface boundary condition
+ `options.bottomBoundary`  bottom boundary condition
+ `options.parameters`  named coefficient parameters

## Returns
+ `evp`  canonical EVP descriptor

## Discussion
