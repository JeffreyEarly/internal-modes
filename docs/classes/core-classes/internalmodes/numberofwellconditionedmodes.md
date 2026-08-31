---
layout: default
title: NumberOfWellConditionedModes
parent: InternalModes
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  NumberOfWellConditionedModes

Estimate how many columns of a mode matrix remain numerically well conditioned.


---

## Declaration
```matlab
 N = NumberOfWellConditionedModes(G,varargin)
```
## Parameters
+ `G`  mode matrix whose columns are ordered by mode number
+ `varargin`  reserved for future options

## Returns
+ `N`  estimated count of usable leading modes

## Discussion

            This function can become a rate limiting step if used for
  each set of returned modes. So a good algorithm is necessary.
  Otherwise we'd just use,
    kappa = InternalModes.ConditionNumberAsFunctionOfModeNumber(G);
    N = find(diff(diff(log10(kappa))) > 1e-2,1,'first')+2;
