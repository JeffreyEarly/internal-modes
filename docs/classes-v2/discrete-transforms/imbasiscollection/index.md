---
layout: default
title: IMBasisCollection
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 1
---

#  IMBasisCollection

Preserve continuous bases and their requested wavenumber-page mapping.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBasisCollection</code></pre></div></div>

## Overview

A collection keeps the original scientific columns and normalization.
It performs no solve, truncation, acceptance decision, or approximate
sharing. `basisIndex` maps requested pages to stored basis objects;
`sourcePage` selects a page within a multi-wavenumber boundary basis.
Supported bases are value objects, so subsequent changes to the input
basis normalization do not change this collection.

```matlab
collection = IMBasisCollection({basis1,basis2},kappa=[k2 k1 k2],basisIndex=[2 1 2]);
G = collection.evaluate(z,variable="G",pages=[3 1]);
```




## Topics
+ Create collections
  + [`IMBasisCollection`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/imbasiscollection.html) Construct a collection without solving or resampling modes.
+ Inspect collections
  + [`bases`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/bases.html) Continuous value bases, in stored-basis order.
  + [`basisIndex`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/basisindex.html) Requested-page to stored-basis mapping.
  + [`kappa`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/kappa.html) Requested horizontal wavenumbers, in requested-page order.
  + [`metadata`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/metadata.html) Scientific column identity and evaluation capabilities per basis.
  + [`sourcePage`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/sourcepage.html) Page within each mapped basis, normally one for scalar solves.
+ Evaluate collections
  + [`evaluate`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/evaluate.html) continuous structures with bounded intermediate allocations.
  + [`evaluateChunks`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/evaluatechunks.html) Stream continuous basis values in bounded sample and page chunks.
+ Assess collections
  + [`assess`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/assess.html) Measure a fixed requested page through the common result contract.
+ Construct projections
  + [`projection`](/internal-modes/classes-v2/discrete-transforms/imbasiscollection/projection.html) Construct a fixed projection for one requested page and variable.


---