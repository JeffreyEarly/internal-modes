---
layout: default
title: leakage
parent: IMProjection
grand_parent: Discrete transforms
nav_order: 12
mathjax: true
---

#  leakage

Measure projected check profiles relative to their positive norms.


---

## Declaration
```matlab
 errors = leakage(self,checkValues,checkNormSquared)
```
## Parameters
+ `checkValues`  finite real or complex sample-by-check matrix
+ `checkNormSquared`  positive squared reference norms, or exact zeros

## Returns
+ `errors`  row of relative projected norms

## Discussion

  Positive norms are supplied explicitly; no check modes are solved here.
  A zero-norm profile reports zero only when its sampled profile is zero;
  otherwise it reports infinity. Values are returned independently per input.
