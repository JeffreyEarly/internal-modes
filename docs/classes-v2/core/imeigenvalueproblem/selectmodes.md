---
layout: default
title: selectModes
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 22
mathjax: true
---

#  selectModes

Select and label retained eigenmodes.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 selection = selectModes(evp,eigenvalues,nModes,context)
```
## Parameters
+ `eigenvalues`  candidate generalized-EVP eigenvalues
+ `nModes`  number of modes to retain
+ `context`  solver or analytical context

## Returns
+ `selection`  selected-mode metadata

## Discussion

  The index policy first classifies candidate eigenvalues using
  boundary metadata, expected null modes, and positive interior
  modes. Retained modes are ordered as boundary-index modes, true
  null modes, then positive interior modes. Their labels define
  the `modeNumber` metadata carried by the resulting `IMBasisSet`.
