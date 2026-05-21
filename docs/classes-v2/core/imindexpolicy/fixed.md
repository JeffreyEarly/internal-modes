---
layout: default
title: fixed
parent: IMIndexPolicy
grand_parent: Classes
nav_order: 9
mathjax: true
---

#  fixed

Create a policy with fixed expected negative and zero counts.


---

## Declaration
```matlab
 policy = IMIndexPolicy.fixed(options)
```
## Parameters
+ `options.expectedNegativeCount`  expected negative count
+ `options.expectedZeroCount`  expected zero count
+ `options.validationMode`  `"error"`, `"warning"`, or `"none"`

## Returns
+ `policy`  initialized index policy

## Discussion
