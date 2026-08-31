---
layout: default
title: transformStateBack
parent: IMGeostrophicTransform
grand_parent: Discrete transforms
nav_order: 10
mathjax: true
---

#  transformStateBack

Reconstruct sampled APV and active endpoint anomalies.


---

## Declaration
```matlab
 [APV,endpointAnomalies] = transformStateBack(transform,options)
```
## Parameters
+ `options.APVCoefficients`  APV coefficient pages
+ `options.zeroAPVCoefficients`  zero-APV coefficient pages
+ `options.zeroAPVCoordinates`  input zero-APV coordinates

## Returns
+ `APV`  sampled APV pages
+ `endpointAnomalies`  endpoint-anomaly pages

## Discussion

  $$
  q=A_{\mathrm i}^{F}\mathbf A_q,\qquad
  \mathbf b^\kappa=\mathsf R_q^\kappa\mathbf A_q^\kappa
  -\frac{f_0}{g\kappa^2}\mathbf A_0^\kappa.
  $$
