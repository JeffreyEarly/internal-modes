---
layout: default
title: rawUz
parent: IMBasisSet
grand_parent: Core
nav_order: 20
mathjax: true
---

#  rawUz

Evaluate raw solved scalar derivatives.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 values = rawUz(basisSet,z)
```
## Parameters
+ `z`  physical coordinate

## Returns
+ `values`  unnormalized scalar derivative values

## Discussion

  This developer utility evaluates the unnormalized physical
  $$z$$-derivative of the solved native modes. The public `uz`
  method applies the active normalization rule after this
  step:
  $$\frac{\partial u_j}{\partial z}(z)=
  \frac{\partial u_j^{\mathrm{raw}}}{\partial z}(z)/s_j,$$
  where $$s_j$$ comes from `normalizationFactors`.
