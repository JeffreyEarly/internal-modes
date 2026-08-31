---
layout: default
title: nativeModes
parent: IMBasisSet
grand_parent: Core
nav_order: 15
mathjax: true
---

#  nativeModes

Native mode columns.

> Developer documentation: this item describes internal implementation details.


---

## Discussion

  These are the columns returned by the numerical solver before
  interpolation onto physical coordinates and before modal
  normalization. Most workflows should call `u(z)` or `uz(z)`
  instead.
