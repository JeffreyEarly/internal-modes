---
layout: default
title: IMProductPlan
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 6
---

#  IMProductPlan

Reserve individual product identities before invoking expensive callbacks.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMProductPlan</code></pre></div></div>

## Overview

The stored inventory is immutable. Execution caches are local to assess;
a later execution cannot reuse stale evaluated fields or reference data.

```matlab
plan = inventory.fixedPlan(prefixCounts=1:8,productBudget=10000);
result = plan.assess(grids,chunkSize=128);
```




## Topics
+ Inspect reservations
  + [`IMProductPlan`](/internal-modes/classes-v2/discrete-transforms/improductplan/improductplan.html) Validate and reserve a complete metadata-only execution plan.
  + [`firstCounts`](/internal-modes/classes-v2/discrete-transforms/improductplan/firstcounts.html) Earliest retained count at which each selected product is examined.
  + [`inventory`](/internal-modes/classes-v2/discrete-transforms/improductplan/inventory.html) Original metadata-only inventory.
  + [`pairs`](/internal-modes/classes-v2/discrete-transforms/improductplan/pairs.html) Two-by-product input column indices for each selected row.
  + [`planningSeconds`](/internal-modes/classes-v2/discrete-transforms/improductplan/planningseconds.html) Inventory construction/planning elapsed seconds.
  + [`prefixCounts`](/internal-modes/classes-v2/discrete-transforms/improductplan/prefixcounts.html) Examined retained-count values, unchanged by assessment decisions.
  + [`productBudget`](/internal-modes/classes-v2/discrete-transforms/improductplan/productbudget.html) Explicit maximum authorized product count.
  + [`reservedProducts`](/internal-modes/classes-v2/discrete-transforms/improductplan/reservedproducts.html) Exact reserved product count, including future structural zeros.
  + [`rows`](/internal-modes/classes-v2/discrete-transforms/improductplan/rows.html) Product-family rows selected in original inventory order.
  + [`selection`](/internal-modes/classes-v2/discrete-transforms/improductplan/selection.html) rule, fixed sparse or bounded all-products validation.
+ Execute assessments
  + [`assess`](/internal-modes/classes-v2/discrete-transforms/improductplan/assess.html) Evaluate only reserved products and reuse preparation across their prefixes.


---