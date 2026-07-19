---
layout: default
title: geometricTransform
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 9
mathjax: true
---

#  geometricTransform

Reference transform using geometric control-volume weights.


---

## Discussion

This transform uses the same points, modes, and normalization as
`transform`, but its weights come only from the physical control
volumes around the fixed points. It provides a baseline for Gram
and transform-quality comparisons.
