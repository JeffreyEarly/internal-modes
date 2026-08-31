---
layout: default
title: IMGeostrophicTransform
parent: IMGeostrophicTransform
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMGeostrophicTransform

Create an APV and zero-APV composition transform.


---

## Declaration
```matlab
 transform = IMGeostrophicTransform(options)
```
## Parameters
+ `options.apvTransform`  generalized-energy APV transform
+ `options.zeroAPVModes`  canonical zero-APV basis
+ `options.g0`  surface acceleration
+ `options.gd`  bottom acceleration
+ `options.muTolerance`  relative singularity tolerance

## Returns
+ `self`  geostrophic composition transform

## Discussion

  Finite endpoint accelerations, including zero, activate the
  corresponding endpoint. Positive infinity makes it inactive.
  Construction rejects modes satisfying

  $$
  \frac{|\mu_\kappa^j|}
  {\kappa^2+|f_0^2/(g h_j)|}\leq\texttt{muTolerance}.
  $$
