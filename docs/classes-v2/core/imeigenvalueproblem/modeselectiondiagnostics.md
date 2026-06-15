---
layout: default
title: modeSelectionDiagnostics
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 10
mathjax: true
---

#  modeSelectionDiagnostics

Summarize negative and zero mode selection.


---

## Declaration
```matlab
 diagnostics = modeSelectionDiagnostics(evp,solver,A)
```
## Parameters
+ `solver`  canonical EVP solver
+ `A`  assembled left matrix

## Returns
+ `diagnostics`  struct with negative and zero mode selection fields

## Discussion

  Negative modes are bounded by `negativeEigenvalueBounds`.
  Zero modes are inferred from the assembled left matrix `A`:
  `zeroModeStatus` is `"present"` when the smallest singular
  value is zero to the rank-style tolerance
  $$\sigma_{\min}(A)\le \max(\mathrm{size}(A))\epsilon(\sigma_{\max}(A)),$$
  `"absent"` when it is larger, and `"unchecked"` when `A` is
  omitted. Mode labels are ordered as
  $$-1,-2,\ldots,\quad 0,\quad 1,2,\ldots.$$
  The returned struct includes `assessmentLevel`,
  `negativeEndpointWeightCount`, `minNegativeEigenvalueCount`,
  `maxNegativeEigenvalueCount`, `zeroModeStatus`,
  `zeroModeCount`, `zeroModeSingularValue`,
  `zeroModeTolerance`, and `reason`.
