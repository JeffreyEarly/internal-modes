---
layout: default
title: selectModes
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 18
mathjax: true
---

#  selectModes

Select and label retained finite-real eigenmodes.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 selection = selectModes(evp,eigenvalues,nModes,solver,A)
```
## Parameters
+ `eigenvalues`  finite real candidate eigenvalues
+ `nModes`  number of retained modes
+ `solver`  canonical solver
+ `A`  assembled left matrix

## Returns
+ `selection`  selected indices and mode numbers

## Discussion

  Mode-selection diagnostics decide when raw negative discrete
  eigenvalues should be retained and whether a zero mode should
  be included. Retained modes are labeled in the order
  $$-1,-2,\ldots,\quad 0,\quad 1,2,\ldots.$$
  The full diagnostics struct is stored in
  `selection.modeSelectionDiagnostics`.
