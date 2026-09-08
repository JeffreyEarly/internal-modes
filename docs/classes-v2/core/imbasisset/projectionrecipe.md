---
layout: default
title: projectionRecipe
parent: IMBasisSet
grand_parent: Core
nav_order: 25
mathjax: true
---

#  projectionRecipe

Bind the scalar signed pairing and continuous evaluation to a recipe.


---

## Parameters
+ `options.variable`  solved scalar variable, "u"

## Returns
+ `recipe`  continuous basis-owned projection recipe

## Discussion

  Scalar targets retain the existing diagonal convention. A signed scalar
  pairing supports Gram measurements but not leakage or quadratic policies.
  Reference integration does not independently qualify the eigenproblem.
