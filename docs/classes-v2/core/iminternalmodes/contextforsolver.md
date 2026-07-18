---
layout: default
title: contextForSolver
parent: IMInternalModes
grand_parent: Core
nav_order: 5
mathjax: true
---

#  contextForSolver

Return the internal-mode coefficient context.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 context = contextForSolver(evp,solver)
```
## Parameters
+ `solver`  canonical solver

## Returns
+ `context`  coefficient context

## Discussion

The returned context extends the canonical context with `N2`,
`f0`, `g`, and `formulation`.
