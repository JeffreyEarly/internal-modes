---
layout: default
title: IMBasisAssessment
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 4
---

#  IMBasisAssessment

Measure a fixed projection without selecting or replacing its columns.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBasisAssessment</code></pre></div></div>

## Overview

Gram, leakage, and supplied-product errors share one result vocabulary.
`measurements` contains one row per requested prefix and quantity.
Measurement status is measured, inconclusive, or notRequested; an
unqualified product reference is inconclusive even if its error is small.
Apply explicit tolerances with `applyPolicy` to obtain a separate decision.

Product recipes supply signed continuous pairings, not full-band
coefficients: each prefix solves its own continuous target Gram system.
The sampled projection keeps its signed metric. Positive coefficient
error norms use the projection's separately supplied majorant.

```matlab
assessment = IMBasisAssessment(projection,prefixCounts=1:8);
decision = assessment.applyPolicy(gramTolerance=1e-2);
```




## Topics
+ Measure projections
  + [`IMBasisAssessment`](/internal-modes/classes-v2/discrete-transforms/imbasisassessment/imbasisassessment.html) Measure Gram quality and optional prepared leakage and products.
  + [`columnKind`](/internal-modes/classes-v2/discrete-transforms/imbasisassessment/columnkind.html) Mode columns or endpoint coordinates; only modes have prefixes.
  + [`costs`](/internal-modes/classes-v2/discrete-transforms/imbasisassessment/costs.html) Measured assessment time, separate from caller construction costs.
  + [`coverage`](/internal-modes/classes-v2/discrete-transforms/imbasisassessment/coverage.html) Actual sampled coverage and omissions, without exhaustive guarantees.
  + [`identity`](/internal-modes/classes-v2/discrete-transforms/imbasisassessment/identity.html) Scientific family, variable, page, and other caller-supplied identity.
  + [`measurements`](/internal-modes/classes-v2/discrete-transforms/imbasisassessment/measurements.html) Quantity, value, status, limiting inputs, and reference provenance.
  + [`prefixCounts`](/internal-modes/classes-v2/discrete-transforms/imbasisassessment/prefixcounts.html) Explicitly examined column counts, including the full requested band.
  + [`projection`](/internal-modes/classes-v2/discrete-transforms/imbasisassessment/projection.html) Original fixed projection, with every requested column preserved.
+ Apply acceptance policies
  + [`applyPolicy`](/internal-modes/classes-v2/discrete-transforms/imbasisassessment/applypolicy.html) Apply explicit tolerances without changing the assessed basis.


---