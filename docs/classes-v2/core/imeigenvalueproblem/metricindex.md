---
layout: default
title: metricIndex
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 13
mathjax: true
---

#  metricIndex

Count negative endpoint directions in the metric.


---

## Declaration
```matlab
 value = metricIndex(evp,options)
```
## Parameters
+ `options.tolerance`  scalar sign tolerance

## Returns
+ `value`  number of active endpoint terms with negative metric weight

## Discussion

  `metricIndex` is the number of active endpoint metric weights
  with negative coefficient, after the supplied sign tolerance.
  It is the finite endpoint part of the indefinite metric index
  used by `negativeEigenvalueBounds`.
