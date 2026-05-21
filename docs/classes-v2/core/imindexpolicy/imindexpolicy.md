---
layout: default
title: IMIndexPolicy
parent: IMIndexPolicy
grand_parent: Classes
nav_order: 1
mathjax: true
---

#  IMIndexPolicy

Create an index policy.


---

## Declaration
```matlab
 policy = IMIndexPolicy(options)
```
## Parameters
+ `options.expectedNegativeCount`  expected negative count
+ `options.expectedZeroCount`  expected zero count
+ `options.indexTolerance`  tolerance for zero classification
+ `options.validationMode`  `"error"`, `"warning"`, or `"none"`
+ `options.boundaryModes`  endpoint boundary-mode slots

## Returns
+ `policy`  initialized index policy

## Discussion
