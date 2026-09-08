---
layout: default
title: IMBasisCollection
parent: IMBasisCollection
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMBasisCollection

Construct a collection without solving or resampling modes.


---

## Parameters
+ `bases`  nonempty cell row of continuous value bases
+ `options.kappa`  requested nonnegative wavenumber row
+ `options.basisIndex`  mapping to stored bases
+ `options.sourcePage`  page within each mapped basis

## Returns
+ `self`  immutable collection of continuous bases

## Discussion

  Multi-page boundary bases require explicit `sourcePage` when
  mapping multiple requested pages. When all mapped bases have
  known wavenumbers, omitted `kappa` is inferred from them.
