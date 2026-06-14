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

  For an active endpoint condition
  $$-[a_i u-b_i(pu')]=\lambda[c_i u-d_i(pu')],$$
  Yassin's indexing uses `z_1` for the bottom and `z_2` for the
  surface, so
  $$D_i=(-1)^{i+1}(a_i d_i-b_i c_i).$$
  Each returned struct has fields `location`, `coefficient`,
  `c`, and `d`, representing the endpoint metric contribution
  $$D_i^{-1}(c_i u-d_i p u_z)^2.$$
