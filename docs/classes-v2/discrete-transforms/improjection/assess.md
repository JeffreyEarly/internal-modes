---
layout: default
title: assess
parent: IMProjection
grand_parent: Discrete transforms
nav_order: 3
mathjax: true
---

#  assess

Measure fixed projection columns without changing them or applying tolerances.


---

## Parameters
+ `options.prefixColumnCounts`  increasing column counts ending at the full band
+ `options.identity`  scientific family and page identity
+ `options.columnKind`  mode or endpoint; endpoints are never prefixed
+ `options.leakage`  sampledValues, normSquared, columnLabels
+ `options.products`  sampledValues, referencePairings, normSquared, inputColumns (2 by n), inputLabels (2 by n), referenceStatus; zero inputColumns denote another family
+ `options.coverage`  caller supplied sampled coverage
+ `options.constructionSeconds`  separate construction cost

## Returns
+ `assessment`  measurements and evidence; checkBasisAssessment applies tolerances

## Discussion

  Return a struct containing projection, identity, columnKind,
  prefixColumnCounts, measurements, coverage, and costs. Measurement rows
  record Gram, leakage, and quadraticAliasing values and reference status.
  Products supply signed referencePairings; every prefix solves its own target
  system. Error magnitudes use the separately supplied positive majorant.
