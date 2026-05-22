---
layout: default
title: active
parent: IMBoundary
grand_parent: Core
nav_order: 1
mathjax: true
---

#  active

Create an active metadata-only boundary condition.


---

## Declaration
```matlab
 boundary = IMBoundary.active(options)
```
## Parameters
+ `options.location`  `"surface"` or `"bottom"`
+ `options.variable`  active variable
+ `options.indexSign`  active-boundary sign
+ `options.indexRank`  number of active directions
+ `options.innerProductTerms`  boundary inner-product terms

## Returns
+ `boundary`  initialized active boundary condition

## Discussion

  Active conditions are already placed because their trace terms
  refer to an endpoint of a partial-depth interval. They declare
  endpoint boundary-mode numbers using `-1` at the surface and
  `-2` at the bottom.
