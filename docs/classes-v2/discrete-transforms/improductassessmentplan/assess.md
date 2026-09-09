---
layout: default
title: assess
parent: IMProductAssessmentPlan
grand_parent: Discrete transforms
nav_order: 2
mathjax: true
---

#  assess

Evaluate only reserved products and reuse preparation across their prefixes.


---

## Parameters
+ `grids`  named sample and explicit reference grids
+ `options.chunkSize`  maximum factor/product columns per temporary batch
+ `options.minimumReciprocalCondition`  explicit coefficient-system guard, default zero

## Returns
+ `result`  measurements, individual evidence, coverage, costs, and separate policies

## Discussion

  `grids` is a struct row with distinct id and finite coordinate column z;
  its first id must be sample. Output prepare(grids) returns projection and
  a cell row of references. Each reference supplies gridId, pairingMatrix,
  normMatrix, targetGramMatrix, majorantGramMatrix, role (primary,
  quadrature, independent), status, and explicit provenance. Exactly one
  primary reference defines errors; other references measure coefficient and
  normalization changes without replacing the assessed scientific modes.
