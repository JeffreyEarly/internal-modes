---
layout: default
title: evaluate
parent: IMOperator
grand_parent: Core
nav_order: 3
mathjax: true
---

#  evaluate

an operator applied to native mode columns.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 values = evaluate(op,solver,nativeModes,z,options)
```
## Parameters
+ `solver`  coordinate-aware internal-mode solver
+ `nativeModes`  native mode columns
+ `z`  physical-coordinate evaluation points
+ `options.context`  framework coefficient context

## Returns
+ `values`  evaluated operator values

## Discussion
