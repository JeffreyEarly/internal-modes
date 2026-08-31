---
layout: default
title: transformStateForward
parent: IMGeostrophicTransform
grand_parent: Discrete transforms
nav_order: 11
mathjax: true
---

#  transformStateForward

Transform an admissible APV and endpoint-anomaly state.


---

## Declaration
```matlab
 [APVCoefficients,zeroAPVCoefficients] = transformStateForward(transform,options)
```
## Parameters
+ `options.APV`  sampled APV pages
+ `options.endpointAnomalies`  endpoint-anomaly pages
+ `options.zeroAPVCoordinates`  output zero-APV coordinates

## Returns
+ `APVCoefficients`  APV coefficient pages
+ `zeroAPVCoefficients`  zero-APV coefficient pages

## Discussion

  $$
  \mathbf A_0^\kappa=-\frac{g\kappa^2}{f_0}
  (\mathbf b^\kappa-\mathsf R_q^\kappa\mathbf A_q^\kappa).
  $$
