---
layout: default
title: IMProductAssessmentPlan
parent: IMProductAssessmentPlan
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMProductAssessmentPlan

Validate declared products and reserve all work before callbacks.


---

## Parameters
+ `factors`  cell row of deferred continuous factors
+ `products`  ordered valid interaction/channel table
+ `outputs`  cell row of deferred physical output projections
+ `options.interactionIds`  supplied interaction IDs; default all
+ `options.retainedCounts`  strictly increasing retained ordinals to assess
+ `options.productBudget`  maximum reservation, including exact zeros
+ `options.selection`  fixedSparse stresses or bounded allProducts control

## Discussion

  Factors declare id, family, labels, ordinals, frequencySigns,
  countRole, evaluate(z,columns,referenceId), and provenance.
  Outputs declare id, family, labels, nondecreasing ordinals,
  countRole, prepare(grids), and provenance. The product table
  requires interactionId, channel, factorA, factorB, and output;
  an optional complex coefficient multiplies each entire product.
