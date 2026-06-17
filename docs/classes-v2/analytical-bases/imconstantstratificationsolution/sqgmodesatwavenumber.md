---
layout: default
title: sqgModesAtWavenumber
parent: IMConstantStratificationSolution
grand_parent: Analytical bases
nav_order: 7
mathjax: true
---

#  sqgModesAtWavenumber

Create exact SQG boundary modes at fixed wavenumber.


---

## Declaration
```matlab
 sqg = sqgModesAtWavenumber(solution,k,options)
```
## Parameters
+ `k`  horizontal wavenumbers
+ `options.boundary`  `"surface"` or `"bottom"`
+ `options.metadata`  additional metadata

## Returns
+ `sqg`  exact SQG basis

## Discussion
