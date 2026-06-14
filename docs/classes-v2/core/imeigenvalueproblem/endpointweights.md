---
layout: default
title: endpointWeights
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 8
mathjax: true
---

#  endpointWeights

Return endpoint metric terms implied by active conditions.


---

## Declaration
```matlab
 weights = endpointWeights(evp,location)
```
## Parameters
+ `location`  `"surface"`, `"bottom"`, or omitted for both endpoints

## Returns
+ `weights`  endpoint metric terms

## Discussion

  Each returned struct has fields `location`, `coefficient`,
  `c`, and `d`, representing
  `coefficient*(c*u-d*p*u_z)^2`.
