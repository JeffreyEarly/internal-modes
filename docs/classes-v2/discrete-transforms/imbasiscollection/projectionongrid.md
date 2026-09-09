---
layout: default
title: projectionOnGrid
parent: IMBasisCollection
grand_parent: Discrete transforms
nav_order: 8
mathjax: true
---

#  projectionOnGrid

Construct a projection on the supplied fixed grid and quadrature.


---

## Parameters
+ `z`  physical sample column
+ `weights`  quadrature weights aligned with z
+ `options.page`  requested page index
+ `options.variable`  sampled variable, formulation by default
+ `options.columns`  explicit array columns, all by default

## Returns
+ `projection`  fixed numerical projection

## Discussion

  Preserve explicit columns, physical signed pairing, and reference provenance.
  No quadrature fitting, acceptance policy, or mode selection is performed.
