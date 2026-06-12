---
layout: default
title: fromBoundarySigns
parent: IMIndexPolicy
grand_parent: Core
nav_order: 10
mathjax: true
---

#  fromBoundarySigns

Create an index policy from active-boundary signs.


---

## Declaration
```matlab
 policy = IMIndexPolicy.fromBoundarySigns(signs,options)
```
## Parameters
+ `signs`  active-boundary signs
+ `options.expectedZeroCount`  expected zero count
+ `options.validationMode`  `"error"`, `"warning"`, or `"none"`

## Returns
+ `policy`  initialized index policy

## Discussion

  Negative signs add one negative-index direction each.
