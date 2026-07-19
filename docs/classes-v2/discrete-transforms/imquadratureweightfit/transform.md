---
layout: default
title: transform
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 18
mathjax: true
---

#  transform

Discrete transform constructed with the fitted weights.


---

## Discussion

This is the production transform associated with `weights`. Its
sample points, retained modes, normalization, metric, and transform
matrices match the quadrature-weight fit.

```matlab
coefficients = weightFit.transform.transformForward(values);
```
