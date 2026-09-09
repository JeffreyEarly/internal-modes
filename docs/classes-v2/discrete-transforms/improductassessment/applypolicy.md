---
layout: default
title: applyPolicy
parent: IMProductAssessment
grand_parent: Discrete transforms
nav_order: 2
mathjax: true
---

#  applyPolicy

Apply explicit product/reference tolerances without reducing counts.


---

## Parameters
+ `options.quadraticTolerance`  nonnegative product-error tolerance
+ `options.referenceTolerance`  maximum qualified reference discrepancy
+ `options.requestedCount`  optional strict examined retained count

## Discussion

  Reference failures remain inconclusive. Qualified product errors
  exceeding tolerance reject. Strict requests refer to an examined
  count and retain all scientific columns supplied by the caller.
