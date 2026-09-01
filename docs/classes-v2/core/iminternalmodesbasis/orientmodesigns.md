---
layout: default
title: orientModeSigns
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 19
mathjax: true
---

#  orientModeSigns

Orient modes so `G` is positive immediately below the surface.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 basisSet = orientModeSigns(basisSet)
```
## Returns
+ `basisSet`  basis set with oriented native mode signs

## Discussion

  A resolved nonzero surface `G` value sets the sign directly.
  When `G` vanishes at the surface, the sign is chosen from
  $$-G_z(z_s)$$, the leading one-sided Taylor coefficient into
  the ocean. The rigid-lid barotropic mode has `G` identically
  zero, so that known `F`-form zero mode uses `F` as a fallback.
  The same sign flip is applied to the coupled `F`/`G` pair.
  `IMInternalModesBasis` is a value class, so callers must keep
  the returned basis set:

  ```matlab
  basisSet = basisSet.orientModeSigns();
  ```
