---
layout: default
title: prefix
parent: IMProjection
grand_parent: Discrete transforms
nav_order: 14
mathjax: true
---

#  prefix

Construct the projection of a leading set of array columns.


---

## Declaration
```matlab
 projection = prefix(self,count)
```
## Parameters
+ `count`  nonnegative number of leading columns

## Returns
+ `projection`  projection on the same fixed sample pairing

## Discussion

  This operation uses array positions. The caller owns whether a leading
  subset is scientifically meaningful; endpoint labels are not mode counts.
