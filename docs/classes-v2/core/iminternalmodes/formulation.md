---
layout: default
title: formulation
parent: IMInternalModes
grand_parent: Core
nav_order: 7
mathjax: true
---

#  formulation

Solved physical variable, `"F"` or `"G"`.


---

## Discussion

  `formulation` tells the canonical solver which variable is the
  native unknown `u`. The complementary variable is evaluated
  diagnostically by `IMInternalModesBasis`. Coefficient handles can
  read this value as `ctx.formulation`.
