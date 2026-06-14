---
layout: default
title: definitenessInfo
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 7
mathjax: true
---

#  definitenessInfo

Check grid-level signs for the canonical coefficients.


---

## Declaration
```matlab
 info = definitenessInfo(evp,solver)
```
## Parameters
+ `solver`  canonical EVP solver

## Returns
+ `info`  struct with sign, metric, and endpoint checks

## Discussion

  This diagnostic certifies the assembled finite-dimensional
  problem on the solver grid. It does not claim continuum signs
  between grid points. The returned struct includes sampled
  minima `pMin`, `qMin`, `rMin`; sign flags `pPositive`,
  `qNonnegative`, and `rPositive`; metric fields
  `endpointWeights`, `metricIndex`, `metricPositive`, and
  `hasDegenerateEndpointMetric`; numerator fields
  `endpointNumeratorNegativeDirections`,
  `endpointNumeratorNonnegative`, `interiorNonnegative`, and
  `qNonnegativeCertified`; and status fields
  `certificationLevel` and `reason`.
