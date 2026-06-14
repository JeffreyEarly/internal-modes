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
+ `options.surfaceBoundary`  surface endpoint condition
+ `options.bottomBoundary`  bottom endpoint condition
+ `options.hasZeroMode`  whether one zero mode should be retained
+ `options.defaultNormalization`  natural normalization
+ `options.normalizations`  named normalization handles
+ `options.parameters`  named coefficient parameters

## Returns
+ `evp`  canonical EVP descriptor

## Discussion
