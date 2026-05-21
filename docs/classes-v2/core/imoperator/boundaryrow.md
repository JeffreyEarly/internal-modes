---
layout: default
title: boundaryRow
parent: IMOperator
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  boundaryRow

Assemble a boundary functional row.


---

## Declaration
```matlab
 row = boundaryRow(op,solver,location,options)
```
## Parameters
+ `solver`  coordinate-aware internal-mode solver
+ `location`  boundary location, `"surface"` or `"bottom"`
+ `options.context`  framework coefficient context

## Returns
+ `row`  assembled boundary row

## Discussion
