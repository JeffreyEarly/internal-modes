---
layout: default
title: orientModeSigns
parent: IMGeostrophicModesBasis
grand_parent: Core
nav_order: 4
mathjax: true
---

#  orientModeSigns

Apply one deterministic sign to each aligned `F`/`G` pair.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 basisSet = orientModeSigns(basisSet)
```
## Returns
+ `basisSet`  basis set with aligned signs oriented
+ `signs`  applied sign for each pair

## Discussion

Resolved surface traces retain the shallow-interior convention.
A trapped signed mode with unresolved surface traces uses its
dominant stored G value; metadata records this phase fallback.
