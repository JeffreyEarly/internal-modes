---
layout: default
title: orientModeSigns
parent: IMBasisSet
grand_parent: Core
nav_order: 22
mathjax: true
---

#  orientModeSigns

Orient scalar modes with a deterministic sign convention.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 basisSet = orientModeSigns(basisSet)
```
## Returns
+ `basisSet`  basis set with oriented native mode signs

## Discussion

  This developer utility flips native mode columns so each
  scalar mode has a deterministic sign based on its largest
  amplitude on the solver integration grid. `IMBasisSet` is a
  value class, so callers must keep the returned basis set:

  ```matlab
  basisSet = basisSet.orientModeSigns();
  ```
