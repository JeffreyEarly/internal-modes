---
layout: default
title: trace
parent: IMBoundary
grand_parent: Core
nav_order: 28
mathjax: true
---

#  trace

Create an endpoint trace descriptor.


---

## Declaration
```matlab
 trace = IMBoundary.trace(variable,options)
```
## Parameters
+ `variable`  variable name
+ `options.derivativeOrder`  physical derivative order

## Returns
+ `trace`  endpoint trace descriptor

## Discussion

  A trace describes which variable value or first derivative is
  evaluated at a boundary endpoint.
