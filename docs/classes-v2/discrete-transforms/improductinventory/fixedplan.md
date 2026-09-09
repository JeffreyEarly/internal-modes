---
layout: default
title: fixedPlan
parent: IMProductInventory
grand_parent: Discrete transforms
nav_order: 5
mathjax: true
---

#  fixedPlan

Reserve a deterministic sparse plan before numerical preparation.


---

## Parameters
+ `options.interactionIds`  caller-selected valid IDs; default all
+ `options.prefixCounts`  increasing requested counts
+ `options.productBudget`  explicit reservation including zeros

## Returns
+ `plan`  metadata-only executable reservation

## Discussion

  Supplied interaction IDs already obey the caller's physical rules. For
  each prefix, retained input factors use low/middle/cutoff ordinal stresses;
  their cumulative union preserves earlier checks. Fixed factors retain every
  column, including every boundary endpoint and frequency-sign coordinate.
