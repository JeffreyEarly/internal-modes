---
layout: default
title: evaluate
parent: IMBasisCollection
grand_parent: Discrete transforms
nav_order: 5
mathjax: true
---

#  evaluate

continuous structures with bounded intermediate allocations.


---

## Parameters
+ `z`  physical coordinate column
+ `options.variable`  scalar variable name
+ `options.derivativeOrder`  supported vertical derivative order
+ `options.pages`  requested-page indices, empty means all
+ `options.columns`  explicit scientific column indices
+ `options.sampleChunkSize`  maximum rows evaluated together
+ `options.pageChunkSize`  maximum pages evaluated together

## Returns
+ `values`  continuous basis values by requested page

## Discussion

  Returns nZ x nColumns x nSelectedPages, preserving page/column order and
  repeated selections. Empty pages/columns selects all. The final output is
  allocated in full; use evaluateChunks to stream without that allocation.
  sampleChunkSize and pageChunkSize bound evaluation temporaries.
