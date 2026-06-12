---
layout: default
title: IMBoundaryWeight
parent: IMBoundaryWeight
grand_parent: Core
nav_order: 1
mathjax: true
---

#  IMBoundaryWeight

Create an endpoint inner-product weight.


---

## Declaration
```matlab
 weight = IMBoundaryWeight(options)
```
## Parameters
+ `options.innerProduct`  inner product receiving the term, `"F"` or `"G"`
+ `options.location`  `"surface"`, `"bottom"`, or `""` for a location-free template
+ `options.coefficient`  scalar or context function handle
+ `options.leftVariable`  left endpoint variable, `"F"` or `"G"`
+ `options.leftDerivativeOrder`  left physical derivative order
+ `options.rightVariable`  right endpoint variable, `"F"` or `"G"`
+ `options.rightDerivativeOrder`  right physical derivative order

## Returns
+ `weight`  initialized boundary weight

## Discussion
