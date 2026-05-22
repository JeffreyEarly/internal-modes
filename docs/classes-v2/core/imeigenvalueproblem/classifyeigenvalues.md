---
layout: default
title: classifyEigenvalues
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 4
mathjax: true
---

#  classifyEigenvalues

Classify eigenvalues using this EVP's index metadata.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 index = classifyEigenvalues(evp,eigenvalues,context)
```
## Parameters
+ `eigenvalues`  generalized-EVP eigenvalues
+ `context`  solver or analytical context

## Returns
+ `index`  index summary structure

## Discussion

  Classification reports the detected boundary-index directions,
  true null modes, positive interior modes, and any validation
  mismatch. Negative boundary directions can be expected by the
  declared boundary policy; they are not automatically treated as
  numerical failures.
