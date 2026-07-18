---
layout: default
title: orientModeSigns
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 12
mathjax: true
---

#  orientModeSigns

Orient modes so the surface `F` value is positive when possible.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 basisSet = orientModeSigns(basisSet)
```
## Returns
+ `basisSet`  basis set with oriented native mode signs

## Discussion

This developer utility applies the internal-mode sign
convention used after numerical solves: prefer a finite
nonzero surface `F` value, fall back to the largest `F`
value on the solver grid, then fall back to `G`. The same
sign flip is applied to the coupled `F`/`G` mode pair.
`IMInternalModesBasis` is a value class, so callers must keep
the returned basis set:

```matlab
basisSet = basisSet.orientModeSigns();
```
