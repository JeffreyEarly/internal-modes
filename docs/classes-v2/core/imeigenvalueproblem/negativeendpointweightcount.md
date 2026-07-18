---
layout: default
title: negativeEndpointWeightCount
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 13
mathjax: true
---

#  negativeEndpointWeightCount

Count negative endpoint metric weights.


---

## Declaration
```matlab
 value = negativeEndpointWeightCount(evp,options)
```
## Parameters
+ `options.tolerance`  scalar sign tolerance

## Returns
+ `value`  number of active endpoint terms with negative metric weight

## Discussion

`negativeEndpointWeightCount` is the number of active
endpoint metric weights with negative coefficient, after the
supplied sign tolerance. The canonical scalar EVP has two
endpoints with at most one active metric weight per endpoint,
so this value is always `0`, `1`, or `2`.
