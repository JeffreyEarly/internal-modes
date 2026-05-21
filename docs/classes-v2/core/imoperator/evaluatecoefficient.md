---
layout: default
title: evaluateCoefficient
parent: IMOperator
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  evaluateCoefficient

Evaluate an operator coefficient on a grid.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 values = IMOperator.evaluateCoefficient(coefficient,z,context)
```
## Parameters
+ `coefficient`  scalar, vector, or function handle coefficient
+ `z`  physical-coordinate evaluation points
+ `context`  framework coefficient context

## Returns
+ `values`  coefficient values matching `z`

## Discussion
