---
layout: default
title: ConditionNumberAsFunctionOfModeNumberForModeIndices
parent: InternalModes
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  ConditionNumberAsFunctionOfModeNumberForModeIndices

Compute mode-matrix condition numbers for selected truncation indices.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [kappa,modeIndices] = ConditionNumberAsFunctionOfModeNumberForModeIndices(G,modeIndices)
```
## Parameters
+ `G`  mode matrix whose columns are ordered by mode number
+ `modeIndices`  truncation indices to evaluate

## Returns
+ `kappa`  condition number at each requested truncation
+ `modeIndices`  echo of the evaluated truncation indices

## Discussion


