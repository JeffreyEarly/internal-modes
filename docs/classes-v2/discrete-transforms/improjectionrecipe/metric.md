---
layout: default
title: metric
parent: IMProjectionRecipe
grand_parent: Discrete transforms
nav_order: 6
mathjax: true
---

#  metric

Construct the signed or majorant pairing on fixed samples.


---

## Parameters
+ `z`  strictly increasing physical sample coordinates
+ `weights`  fixed quadrature weights
+ `options.majorant`  use absolute endpoint coefficients

## Returns
+ `matrix`  nZ by nZ sample metric

## Discussion

  Endpoint traces must be represented by supplied variable
  samples. Derivative and companion-variable traces are explicit
  unsupported capabilities, never inferred from the basis span.
  Quadrature weights are unchanged; negative quadrature weights
  can make the sampled majorant indefinite.
