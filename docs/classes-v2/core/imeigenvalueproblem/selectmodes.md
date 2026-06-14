---
layout: default
title: selectModes
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 22
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

  Certified negative-count bounds decide when raw negative
  discrete eigenvalues should be retained. With a positive metric
  and nonnegative quadratic form, negative discrete eigenvalues
  are ignored during mode selection.
