---
layout: default
title: boundaryConditions
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 3
mathjax: true
---

#  boundaryConditions

Placed boundary-condition array.


---

## Discussion

  The array stores location-aware `IMBoundary` values. Standard
  factories accept location-free `surfaceBoundary` and `bottomBoundary`
  laws and place them on the solved variable. During assembly each
  placed boundary replaces the appropriate boundary row in `A` and `B`
  and contributes endpoint metadata for indexing and normalization.
