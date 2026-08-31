---
layout: default
title: definitenessDiagnostics
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 5
mathjax: true
---

#  definitenessDiagnostics

Assess the norm and energy signs that control negative modes.


---

## Declaration
```matlab
 diagnostics = definitenessDiagnostics(evp,solver)
```
## Parameters
+ `solver`  canonical EVP solver

## Returns
+ `diagnostics`  struct with norm and energy checks

## Discussion

  On the solver grid, the canonical EVP has a norm determined
  by the interior weight `r` and endpoint weights, and an energy
  determined by the interior `p`, `q` terms plus endpoint
  contributions. This diagnostic checks the sampled coefficients
  and endpoint terms that determine whether the norm is positive
  and whether the energy is nonnegative. These are grid-level
  checks for the discretized problem, not continuum guarantees
  between grid points.

  The returned struct includes sampled minima `pMin`, `qMin`,
  `rMin`; sign flags `pPositive`, `qNonnegative`, and
  `rPositive`; norm fields `endpointWeights`,
  `negativeEndpointWeightCount`, `metricPositive`, and
  `hasDegenerateEndpointMetric`; energy fields
  `endpointNumeratorNegativeDirections`,
  `endpointNumeratorNonnegative`, `interiorNonnegative`, and
  `quadraticFormNonnegative`; and status fields
  `assessmentLevel` and `reason`.
