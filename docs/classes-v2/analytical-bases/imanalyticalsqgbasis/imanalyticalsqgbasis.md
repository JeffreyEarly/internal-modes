---
layout: default
title: IMAnalyticalSQGBasis
parent: IMAnalyticalSQGBasis
grand_parent: Analytical bases
nav_order: 1
mathjax: true
---

#  IMAnalyticalSQGBasis

Create an exact SQG basis.


---

## Declaration
```matlab
 sqg = IMAnalyticalSQGBasis(options)
```
## Parameters
+ `options.solution`  analytical solution family
+ `options.k`  horizontal wavenumbers
+ `options.boundary`  active boundary
+ `options.N2`  buoyancy frequency squared function
+ `options.psiFunction`  exact streamfunction evaluator
+ `options.metadata`  creation metadata

## Returns
+ `sqg`  exact SQG basis

## Discussion
