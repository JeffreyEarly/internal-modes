---
layout: default
title: transformBack
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 16
mathjax: true
---

#  transformBack

Transform modal coefficients back to sampled profiles.


---

## Declaration
```matlab
 values = transformBack(transform,coefficients)
```
## Parameters
+ `coefficients`  `nModes`-by-`nProfiles` coefficient array with rows aligned to `modeNumber`

## Returns
+ `values`  `nSamples`-by-`nProfiles` reconstructed profile array sampled on `z`

## Discussion

  For $$n_p$$ coefficient sets arranged as
  $$A\in\mathbb{R}^{n_m\times n_p}$$, this method returns
  $$\widehat{X}\in\mathbb{R}^{n_z\times n_p}$$ using

  $$
  \widehat{X}=A_{\mathrm i}A=\Phi A.
  $$

  Each coefficient column is transformed independently. For
  coefficients obtained from `transformForward`, the result is
  the sampled-space projection
  $$\widehat{X}=A_{\mathrm i}A_{\mathrm f}X$$, not generally the
  original profile unless it lies in the retained modal subspace.

  ```matlab
  values = transform.transformBack(coefficients);
  valuesByMatrix = transform.inverseMatrix*coefficients;
  ```
