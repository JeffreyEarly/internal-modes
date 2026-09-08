---
layout: default
title: applyPolicy
parent: IMBasisAssessment
grand_parent: Discrete transforms
nav_order: 2
mathjax: true
---

#  applyPolicy

Apply explicit tolerances without changing the assessed basis.


---

## Declaration
```matlab
 decision = applyPolicy(assessment,options)
```
## Discussion

  Disabled quantities do not enter the decision. Missing or
  unqualified enabled measurements are inconclusive. A measured
  failure rejects even when another enabled quantity is inconclusive.
  Prefix recommendations are diagnostic only; requestedColumnCount
  and the stored projection are never reduced.
