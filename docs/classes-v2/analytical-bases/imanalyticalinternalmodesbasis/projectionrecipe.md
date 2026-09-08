---
layout: default
title: projectionRecipe
parent: IMAnalyticalInternalModesBasis
grand_parent: Analytical bases
nav_order: 26
mathjax: true
---

#  projectionRecipe

Bind an analytical aligned basis without requiring a numerical solver.


---

## Parameters
+ `options.variable`  "F" or "G", the solved formulation by default

## Returns
+ `recipe`  continuous basis-owned projection recipe

## Discussion

  Gram targets use the existing 1024-point trapezoidal reference rule;
  analytical evaluation does not make this quadrature convergence-qualified.
