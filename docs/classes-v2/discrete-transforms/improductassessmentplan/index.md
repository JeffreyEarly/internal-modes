---
layout: default
title: IMProductAssessmentPlan
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 3
---

#  IMProductAssessmentPlan

Reserve individual product identities before invoking expensive callbacks.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMProductAssessmentPlan</code></pre></div></div>

## Overview

The stored inventory is immutable. Execution caches are local to assess;
a later execution cannot reuse stale evaluated fields or reference data.

```matlab
plan = IMProductAssessmentPlan(factors,products,outputs,retainedCounts=1:8,productBudget=10000);
result = plan.assess(grids,chunkSize=128);
```




## Topics
+ Inspect reservations
  + [`IMProductAssessmentPlan`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/improductassessmentplan.html) Validate declared products and reserve all work before callbacks.
  + [`factors`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/factors.html) Deferred continuous factor metadata and evaluators.
  + [`firstCounts`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/firstcounts.html) Earliest retained count at which each selected product is examined.
  + [`outputs`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/outputs.html) Deferred output projection metadata and preparation callbacks.
  + [`pairs`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/pairs.html) Two-by-product input column indices for each selected row.
  + [`planningSeconds`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/planningseconds.html) Inventory construction/planning elapsed seconds.
  + [`productBudget`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/productbudget.html) Explicit maximum authorized product count.
  + [`products`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/products.html) Ordered declared interaction and channel identities.
  + [`reservedProducts`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/reservedproducts.html) Exact reserved product count, including future structural zeros.
  + [`retainedCounts`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/retainedcounts.html) Examined retained-count values, unchanged by assessment decisions.
  + [`rows`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/rows.html) Product-family rows selected in original inventory order.
  + [`selection`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/selection.html) rule, fixed sparse or bounded all-products validation.
+ Execute assessments
  + [`assess`](/internal-modes/classes-v2/discrete-transforms/improductassessmentplan/assess.html) Evaluate only reserved products and reuse preparation across their prefixes.


---