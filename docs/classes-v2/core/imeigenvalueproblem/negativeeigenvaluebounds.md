---
layout: default
title: negativeEigenvalueBounds
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 14
mathjax: true
---

#  negativeEigenvalueBounds

Bound negative eigenvalues using grid-level certification.


---

## Declaration
```matlab
 bounds = negativeEigenvalueBounds(evp,solver,A)
```
## Parameters
+ `solver`  canonical EVP solver
+ `A`  assembled left matrix, used for the zero-eigenvalue check

## Returns
+ `bounds`  struct with min/max counts and a reason

## Discussion

  The returned counts describe how many negative eigenvalues are
  certified by the discretized canonical problem, rather than by
  raw negative finite-real eigenvalues alone. The returned struct
  includes `certificationLevel`, `negativeEndpointWeightCount`,
  `zeroEigenvalueStatus`, `minNegativeEigenvalueCount`,
  `maxNegativeEigenvalueCount`, and `reason`.
  `maxNegativeEigenvalueCount` may be the string `"unknown"`
  when coefficient signs or endpoint determinants cannot be
  certified on the grid.
