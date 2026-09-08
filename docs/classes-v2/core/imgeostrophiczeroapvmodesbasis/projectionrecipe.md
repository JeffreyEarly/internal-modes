---
layout: default
title: projectionRecipe
parent: IMGeostrophicZeroAPVModesBasis
grand_parent: Core
nav_order: 17
mathjax: true
---

#  projectionRecipe

Report the unavailable scalar F/G projection for boundary-response bases.


---

## Parameters
+ `options.variable`  requested scalar variable, "F" or "G"

## Returns
+ `recipe`  unavailable recipe with an explicit reason

## Discussion

  Endpoint response and energy matrices are coefficient-space forms. They
  do not define projection of arbitrary scalar F or G samples. Endpoint
  labels remain intact; no modal-prefix truncation is inferred.
