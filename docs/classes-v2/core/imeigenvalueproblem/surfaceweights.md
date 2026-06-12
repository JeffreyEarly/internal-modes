---
layout: default
title: surfaceWeights
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 26
mathjax: true
---

#  surfaceWeights

Surface endpoint weights implied by the surface boundary law.


---

## Discussion

  The returned `IMBoundaryWeight` array is derived from
  `surfaceBoundary` and `formulation`. Changing the surface law or
  formulation changes this read-only view automatically.
