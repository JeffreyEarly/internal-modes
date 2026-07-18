---
layout: default
title: project
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 12
mathjax: true
---

#  project

Return retained modal coefficients for sampled profiles.


---

## Declaration
```matlab
 coefficients = project(transform,values)
```
## Parameters
+ `values`  sampled profiles with rows aligned to `z`

## Returns
+ `coefficients`  retained modal coefficients

## Discussion

This method applies `forwardMatrix`:

$$
\mathbf{a}=A_{\mathrm f}\mathbf{x}.
$$
