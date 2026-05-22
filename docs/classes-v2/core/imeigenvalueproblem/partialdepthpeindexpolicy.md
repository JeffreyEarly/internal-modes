---
layout: default
title: partialDepthPEIndexPolicy
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 20
mathjax: true
---

#  partialDepthPEIndexPolicy

Return the partial-depth potential-energy index policy.


---

## Declaration
```matlab
 policy = IMEigenvalueProblem.partialDepthPEIndexPolicy(options)
```
## Parameters
+ `options.boundarySign`  `"positive"` or `"negative"`
+ `options.validationMode`  `"error"`, `"warning"`, or `"none"`

## Returns
+ `policy`  corresponding index policy

## Discussion

  Partial-depth potential-energy forms can declare endpoint
  directions through active boundary metadata. With
  `boundarySign="positive"`, endpoint contributions use the
  positive boundary convention. With `boundarySign="negative"`,
  endpoint contributions use the negative boundary convention and
  the resulting negative index directions are expected by the
  policy rather than treated as numerical failures.
