---
layout: default
title: evaluateChunks
parent: IMBasisCollection
grand_parent: Discrete transforms
nav_order: 6
mathjax: true
---

#  evaluateChunks

Stream continuous basis values in bounded sample and page chunks.


---

## Parameters
+ `z`  physical coordinate column
+ `consumer`  callback receiving values, rows, positions, and page IDs
+ `options.variable`  scalar variable name, first available by default
+ `options.derivativeOrder`  supported vertical derivative order
+ `options.pages`  requested collection pages, empty means all
+ `options.columns`  explicit scientific column indices
+ `options.sampleChunkSize`  maximum coordinate rows per callback
+ `options.pageChunkSize`  maximum requested pages per callback

## Discussion

  Calls `consumer(values,rowIndices,pagePositions,requestedPages)`, where
  pagePositions indexes the supplied `pages` vector and requestedPages holds
  the corresponding collection page IDs. This preserves repeated page IDs
  without ambiguity. Values have shape nRows x nColumns x nChunkPages.
  No full output array is allocated. Validation precedes the first callback.
  Compatible native spectral bases share evaluation matrices exactly;
  analytical bases and custom subclasses retain their own evaluators.
  Custom evaluators control their internal allocations; chunk bounds apply
  to callbacks and the built-in evaluation paths, not arbitrary user code.
