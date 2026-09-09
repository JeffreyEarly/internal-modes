---
layout: default
title: compareCoverage
parent: IMProductAssessment
grand_parent: Discrete transforms
nav_order: 4
mathjax: true
---

#  compareCoverage

Detect failures missed by this plan using a separate bounded control.


---

## Parameters
+ `control`  independently selected or all-products assessment
+ `options.quadraticTolerance`  explicit failure threshold
+ `options.referenceTolerance`  required reference qualification

## Discussion

  A control must use the same declared inventory and counts, with
  qualified references. This reports observed missed products;
  no failure can be detected outside both examined inventories.
