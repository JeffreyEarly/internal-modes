---
layout: default
title: projection
parent: IMProjectionRecipe
grand_parent: Discrete transforms
nav_order: 8
mathjax: true
---

#  projection

Construct a projection for explicit columns on an unchanged grid.


---

## Parameters
+ `z`  physical sample coordinates
+ `weights`  fixed quadrature weights
+ `options.columns`  explicit column selection, all by default

## Returns
+ `projection`  numerical projection and measurements
