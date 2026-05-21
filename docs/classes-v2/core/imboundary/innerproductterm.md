---
layout: default
title: innerProductTerm
parent: IMBoundary
grand_parent: Classes
nav_order: 17
mathjax: true
---

#  innerProductTerm

Create a boundary inner-product trace-pair term.


---

## Declaration
```matlab
 term = IMBoundary.innerProductTerm(innerProductVariable,location,coefficient,leftTrace,rightTrace)
```
## Parameters
+ `innerProductVariable`  variable whose inner product receives the term
+ `location`  boundary location
+ `coefficient`  scalar or context function handle
+ `leftTrace`  trace evaluated for the left mode
+ `rightTrace`  trace evaluated for the right mode

## Returns
+ `term`  boundary inner-product term

## Discussion

  The term contributes `coefficient*leftTrace_i*rightTrace_j`
  at a boundary endpoint to the named variable's inner product.
