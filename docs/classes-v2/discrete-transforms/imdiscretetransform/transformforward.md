---
layout: default
title: transformForward
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 17
mathjax: true
---

#  transformForward

Transform sampled profiles forward to modal coefficients.


---

## Declaration
```matlab
 coefficients = transformForward(transform,values)
```
## Parameters
+ `values`  `nSamples`-by-`nProfiles` sampled profile array with rows aligned to `z`

## Returns
+ `coefficients`  `nModes`-by-`nProfiles` retained modal coefficient array

## Discussion

For $$n_p$$ profiles arranged as
$$X\in\mathbb{R}^{n_z\times n_p}$$, this method returns
$$A\in\mathbb{R}^{n_m\times n_p}$$ using

$$
A=A_{\mathrm f}X
=\left(A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}\right)^{-1}
A_{\mathrm i}^\mathsf{T}WX.
$$

Each input column is transformed independently. The resulting
coefficients minimize the sampled quadratic residual when
$$W$$ is positive definite and, more generally, satisfy

$$
A_{\mathrm i}^\mathsf{T}W
\left(X-A_{\mathrm i}A\right)=0.
$$

```matlab
coefficients = transform.transformForward(values);
coefficientsByMatrix = transform.forwardMatrix*values;
```
