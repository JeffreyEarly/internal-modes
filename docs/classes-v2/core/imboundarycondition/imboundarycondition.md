---
layout: default
title: IMBoundaryCondition
parent: IMBoundaryCondition
grand_parent: Core
nav_order: 1
mathjax: true
---

#  IMBoundaryCondition

Create a scalar endpoint condition.


---

## Declaration
```matlab
 boundary = IMBoundaryCondition(options)
```
## Parameters
+ `options.a`  value coefficient on the left
+ `options.b`  flux coefficient on the left
+ `options.c`  value coefficient on the eigenvalue side
+ `options.d`  flux coefficient on the eigenvalue side

## Returns
+ `boundary`  endpoint condition

## Discussion
