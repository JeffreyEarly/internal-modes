---
layout: default
title: selectModes
parent: IMEigenvalueProblem
grand_parent: Classes
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

  Boundary conditions provide boundary-mode index metadata.
  `nNullModes` provides expected zero-eigenvalue null modes. The
  selected modes are ordered as boundary modes, null modes, then
  positive interior modes.
