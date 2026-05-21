---
layout: default
title: fromBoundaryConditions
parent: IMIndexPolicy
grand_parent: Classes
nav_order: 10
mathjax: true
---

#  fromBoundaryConditions

Create an index policy from boundary-condition index metadata.


---

## Declaration
```matlab
 policy = IMIndexPolicy.fromBoundaryConditions(boundaryConditions,options)
```
## Parameters
+ `boundaryConditions`  boundary-condition array
+ `options.expectedZeroCount`  additional expected zero count
+ `options.validationMode`  `"error"`, `"warning"`, or `"none"`

## Returns
+ `policy`  initialized index policy

## Discussion

  Placed boundary conditions contribute their negative and zero
  index counts. Conditions with unknown compatible inner-product
  terms do not contribute expected counts.
