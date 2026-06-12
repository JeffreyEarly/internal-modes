---
layout: default
title: at
parent: IMBoundaryWeight
grand_parent: Core
nav_order: 2
mathjax: true
---

#  at

Place location-free weights at a physical endpoint.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 weights = at(weights,location)
```
## Parameters
+ `location`  `"surface"` or `"bottom"`

## Returns
+ `weights`  placed boundary weights

## Discussion

  Explicitly located weights must already match the requested
  endpoint. Location-free weights receive the endpoint location;
  bottom placement flips the coefficient sign.
