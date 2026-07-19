---
layout: default
title: geometricWeights
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 10
mathjax: true
---

#  geometricWeights

Geometric control-volume weights aligned with the fixed points.


---

## Discussion

For adjacent points, control-volume edges lie at their midpoints;
the outer edges are the basis-set domain boundaries. Therefore
these weights are geometric widths and sum to the full depth.
They satisfy `geometricWeights = geometricTransform.weights`.
