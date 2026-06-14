---
layout: default
title: metricWeight
parent: IMBoundaryCondition
grand_parent: Core
nav_order: 10
mathjax: true
---

#  metricWeight

Return the endpoint metric weight.


---

## Declaration
```matlab
 value = metricWeight(boundary,location)
```
## Parameters
+ `location`  `"surface"` or `"bottom"`

## Returns
+ `value`  coefficient multiplying `(c*u-d*p*u_z)^2`

## Discussion
