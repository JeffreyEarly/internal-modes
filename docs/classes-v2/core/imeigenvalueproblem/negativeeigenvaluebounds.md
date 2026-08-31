---
layout: default
title: negativeEigenvalueBounds
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 12
mathjax: true
---

#  negativeEigenvalueBounds

Bound negative eigenvalues using a grid-level assessment.


---

## Declaration
```matlab
 bounds = negativeEigenvalueBounds(evp,solver,A)
```
## Parameters
+ `solver`  canonical EVP solver
+ `A`  assembled left matrix, used for the zero-mode check

## Returns
+ `bounds`  struct with min/max counts and a reason

## Discussion

  The returned counts describe how many negative eigenvalues are
  supported by the discretized canonical problem, rather than by
  raw negative finite-real eigenvalues alone. Exact negative
  counts require the zero mode to be absent. The returned struct
  includes `assessmentLevel`, `negativeEndpointWeightCount`,
  `zeroModeStatus`, `minNegativeEigenvalueCount`,
  `maxNegativeEigenvalueCount`, and `reason`.
  `maxNegativeEigenvalueCount` may be the string `"unknown"`
  when coefficient signs or endpoint determinants cannot be
  assessed on the grid.
