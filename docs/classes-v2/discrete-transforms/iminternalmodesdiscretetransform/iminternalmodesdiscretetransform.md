---
layout: default
title: IMInternalModesDiscreteTransform
parent: IMInternalModesDiscreteTransform
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMInternalModesDiscreteTransform

Create an aligned internal-mode transform from prepared data.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 transform = IMInternalModesDiscreteTransform(options)
```
## Discussion

  `channelData.F` and `channelData.G` are scalar structs with
  fields `available`, `reason`, `activeModeMask`,
  `metricMatrix`, `targetGramMatrix`, and
  `targetMajorantGramMatrix`. Diagnostic inverse
  matrices and endpoint traces are supplied for both variables,
  including variables that have no direct sampled projection.
