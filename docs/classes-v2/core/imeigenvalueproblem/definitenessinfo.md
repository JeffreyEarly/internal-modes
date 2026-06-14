---
layout: default
title: definitenessInfo
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 6
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
  between grid points.
