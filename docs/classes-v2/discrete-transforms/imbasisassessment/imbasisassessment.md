---
layout: default
title: IMBasisAssessment
parent: IMBasisAssessment
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMBasisAssessment

Measure Gram quality and optional prepared leakage and products.


---

## Declaration
```matlab
 assessment = IMBasisAssessment(projection,options)
```
## Discussion

  `leakage` contains sampledValues, normSquared, and columnLabels.
  Check columns must belong to this same scientific family and page.
  Labels matching retained columns are excluded at each prefix.
  `products` contains sampledValues, referencePairings,
  normSquared, inputColumns (2-by-nProducts in the selected output
  column order; zero denotes an input
  outside this output family), inputLabels (2-by-nProducts), and
  referenceStatus (qualified or inconclusive/unverified).
  Optional referenceProvenance describes the independent evidence.
