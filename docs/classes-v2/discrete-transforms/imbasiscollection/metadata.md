---
layout: default
title: metadata
parent: IMBasisCollection
grand_parent: Discrete transforms
nav_order: 8
mathjax: true
---

#  metadata

Scientific column identity and evaluation capabilities per basis.


---

## Type
+ Class: `struct`
+ Size: `(1,:)`

## Discussion

  `derivativeOrders` is a cell row aligned with `variables`. These
  describe evaluation only; no projection capability is implied.
  Mode labels, endpoint labels, and wavenumber pages remain separate.
  Frequency signs are not introduced by a vertical basis collection.
